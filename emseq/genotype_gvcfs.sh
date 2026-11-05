#!/bin/bash
#SBATCH --job-name=GenotypeGVCF
#SBATCH --array=0-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=variant-calling/logs/Genotype_%A_%a.out
#SBATCH --error=variant-calling/logs/Genotype_%A_%a.err

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mamba_abner_BC

module load gatk/4.1.8.0

REF="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/GRCm39.fa"

GVCF_DIR="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/variant-calling/gvcf"
VCF_DIR="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/variant-calling/vcf"

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

gatk --java-options "-Xmx20g" GenotypeGVCFs \
    -R "$REF" \
    -V "${GVCF_DIR}/${SAMPLE}.g.vcf.gz" \
    -O "${VCF_DIR}/${SAMPLE}.raw.vcf.gz" \
    --tmp-dir "${SLURM_TMPDIR:-/tmp}"

tabix -f -p vcf "${VCF_DIR}/${SAMPLE}.raw.vcf.gz"

echo "Sample ${SAMPLE}: $(zcat ${VCF_DIR}/${SAMPLE}.raw.vcf.gz | grep -vc '^#') variants" >> "${VCF_DIR}/variant_counts.txt"