#!/bin/bash
#SBATCH --job-name=WGS_align
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --output=wgs_alignment_%j.out
#SBATCH --error=wgs_alignment_%j.err

set -euo pipefail

# Activate conda environment
source ~/.bashrc
conda activate mamba_abner_BC

# Reference
REF=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/GRCm39.fa

# FASTQ directory
FASTQ=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/fastq

# Output directory
OUT=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/wgs-alignment

mkdir -p "$OUT"

echo "Starting WGS alignment"
date

# -----------------------------
# AMPK-KO-CD
# -----------------------------

bwa mem -t 16 \
    -R '@RG\tID:YRLP_0066\tSM:YRLP_0066\tPL:ILLUMINA' \
    "$REF" \
    "$FASTQ/YRLP_0066_1_CO_Whole_C1_NFPFG_A48660_237232LT4_ACAGTGAC_L007_R1_001.fastq.gz" \
    "$FASTQ/YRLP_0066_1_CO_Whole_C1_NFPFG_A48660_237232LT4_ACAGTGAC_L007_R2_001.fastq.gz" \
    | samtools sort -@ 8 -m 4G \
    -o "$OUT/YRLP_0066.sorted.bam" -

samtools index "$OUT/YRLP_0066.sorted.bam"

# -----------------------------
# AMPK-WT-CD
# -----------------------------

bwa mem -t 16 \
    -R '@RG\tID:YRLP_0067\tSM:YRLP_0067\tPL:ILLUMINA' \
    "$REF" \
    "$FASTQ/YRLP_0067_1_CO_Whole_C1_NFPFG_A48661_237232LT4_TGACCAAC_L007_R1_001.fastq.gz" \
    "$FASTQ/YRLP_0067_1_CO_Whole_C1_NFPFG_A48661_237232LT4_TGACCAAC_L007_R2_001.fastq.gz" \
    | samtools sort -@ 8 -m 4G \
    -o "$OUT/YRLP_0067.sorted.bam" -

samtools index "$OUT/YRLP_0067.sorted.bam"

# -----------------------------
# CT26
# -----------------------------

bwa mem -t 16 \
    -R '@RG\tID:YRLP_0068\tSM:YRLP_0068\tPL:ILLUMINA' \
    "$REF" \
    "$FASTQ/YRLP_0068_1_CL_Whole_C1_NFPFG_A48662_237232LT4_CGATGTAC_L007_R1_001.fastq.gz" \
    "$FASTQ/YRLP_0068_1_CL_Whole_C1_NFPFG_A48662_237232LT4_CGATGTAC_L007_R2_001.fastq.gz" \
    | samtools sort -@ 8 -m 4G \
    -o "$OUT/YRLP_0068.sorted.bam" -

samtools index "$OUT/YRLP_0068.sorted.bam"

echo "Alignment complete"
date