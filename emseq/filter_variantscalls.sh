#!/bin/bash
#SBATCH --job-name=VF
#SBATCH --array=0-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --output=logs/VF_%A_%a.out
#SBATCH --error=logs/VF_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -eo pipefail

module load gatk          
set -u

which gatk
gatk --version 2>&1 | head -n 2

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REF=$BASE/genotype-reference/GRCm39.fa
MIN_GQ=20

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
MIN_DPS=(10 12 14)
MAX_DPS=(55 70 80)

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
MIN_DP=${MIN_DPS[$SLURM_ARRAY_TASK_ID]}
MAX_DP=${MAX_DPS[$SLURM_ARRAY_TASK_ID]}

IN=$BASE/variant-calling/vcf/${S}.raw.vcf.gz
OUT=$BASE/variant-calling/filtered
mkdir -p $OUT logs
cd $OUT

# 1. Split SNPs and indels
gatk SelectVariants -R $REF -V $IN --select-type-to-include SNP   -O ${S}.snp.vcf.gz
gatk SelectVariants -R $REF -V $IN --select-type-to-include INDEL -O ${S}.indel.vcf.gz

# 2. Mark bad sites (SNPs)
gatk VariantFiltration -R $REF -V ${S}.snp.vcf.gz \
-filter "QD < 2.0" --filter-name QD2 \
-filter "QUAL < 30.0" --filter-name QUAL30 \
-filter "SOR > 3.0" --filter-name SOR3 \
-filter "FS > 60.0" --filter-name FS60 \
-filter "MQ < 40.0" --filter-name MQ40 \
-filter "MQRankSum < -12.5" --filter-name MQRS \
-filter "ReadPosRankSum < -8.0" --filter-name RPRS \
--genotype-filter-expression "DP < ${MIN_DP} || DP > ${MAX_DP} || GQ < ${MIN_GQ}" --genotype-filter-name lowDPGQ \
-O ${S}.snp.filt.vcf.gz

# 2b. Mark bad sites (indels)
gatk VariantFiltration -R $REF -V ${S}.indel.vcf.gz \
-filter "QD < 2.0" --filter-name QD2 \
-filter "QUAL < 30.0" --filter-name QUAL30 \
-filter "FS > 200.0" --filter-name FS200 \
-filter "ReadPosRankSum < -20.0" --filter-name RPRS \
--genotype-filter-expression "DP < ${MIN_DP} || DP > ${MAX_DP} || GQ < ${MIN_GQ}" --genotype-filter-name lowDPGQ \
-O ${S}.indel.filt.vcf.gz

# 3. Drop flagged sites, blank out low-confidence genotypes
for T in snp indel; do
gatk SelectVariants -R $REF -V ${S}.${T}.filt.vcf.gz \
--exclude-filtered --set-filtered-gt-to-nocall -O ${S}.${T}.pass.vcf.gz
done

# 4. Recombine
gatk MergeVcfs -I ${S}.snp.pass.vcf.gz -I ${S}.indel.pass.vcf.gz -O ${S}.pass.vcf.gz

# Put bcftools/tabix env on the PATH for the remaining steps
export PATH=/home/abportillo/.local/share/mamba/envs/bcftools_env/bin:$PATH
bcftools --version | head -n 1

# 5. Normalize: split multiallelics, left-align indels, remove '*' alleles
bcftools norm -m -any -f $REF ${S}.pass.vcf.gz -Ou \
| bcftools view -e 'ALT="*"' -Oz -o ${S}.pass.norm.vcf.gz
tabix -f -p vcf ${S}.pass.norm.vcf.gz

# 6. Keep only homozygous ALT sites for reference building
bcftools view -i 'GT="1/1"' ${S}.pass.norm.vcf.gz -Oz -o ${S}.homalt.vcf.gz
tabix -f -p vcf ${S}.homalt.vcf.gz

# 7. QC summaries
bcftools stats ${S}.pass.norm.vcf.gz > ${S}.pass.stats.txt
bcftools stats ${S}.homalt.vcf.gz    > ${S}.homalt.stats.txt