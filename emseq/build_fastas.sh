#!/bin/bash
#SBATCH --job-name=sample_refs
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/CONS_%A_%a.out
#SBATCH --error=logs/CONS_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -eo pipefail
export PATH=/home/abportillo/.local/share/mamba/envs/bcftools_env/bin:$PATH
set -u

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REF=$BASE/genotype-reference/GRCm39.fa
VCFDIR=$BASE/variant-calling/filtered
OUT=$BASE/genotype-reference/sample-refs
mkdir -p $OUT logs

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

bcftools consensus -f $REF -H 1 -c $OUT/${S}.chain \
  $VCFDIR/${S}.homalt.vcf.gz > $OUT/${S}.fa 2> $OUT/${S}.consensus.log

# applied-variant count vs VCF record count
echo "VCF records: $(bcftools view -H $VCFDIR/${S}.homalt.vcf.gz | wc -l)"
cat $OUT/${S}.consensus.log