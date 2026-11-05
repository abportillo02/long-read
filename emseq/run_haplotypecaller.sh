#!/bin/bash
#SBATCH --job-name=HC
#SBATCH --array=0-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=variant-calling/logs/HC_%A_%a.out
#SBATCH --error=variant-calling/logs/HC_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=abportillo@coh.org


set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mamba_abner_BC

module load gatk/4.1.8.0

REF="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/GRCm39.fa"

BAM_DIR="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/wgs-alignment"
OUT_DIR="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/variant-calling/gvcf"

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

BAM="${BAM_DIR}/${SAMPLE}.markdup.bam"
GVCF="${OUT_DIR}/${SAMPLE}.g.vcf.gz"

gatk --java-options "-Xmx28g" HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$GVCF" \
    -ERC GVCF \
    --native-pair-hmm-threads 8 \
    --tmp-dir "${SLURM_TMPDIR:-/tmp}"

tabix -p vcf "$GVCF"
