#!/bin/bash
#SBATCH --job-name=kzfp_isoform_calling
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=kzfp_isoform_calling_%j.out
#SBATCH --error=kzfp_isoform_calling_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -uo pipefail

# --- Environment (samtools + minimap2 live here) ---
set +u
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC
set -u

THREADS=${SLURM_CPUS_PER_TASK:-16}

# --- Paths ---
GENOME=/home/abportillo/genomes/hg38/hg38_p14.fa
OUTDIR=/home/abportillo/github_repo/long-read/results/pilot_kzfp
BED=$OUTDIR/panel.genes.bed
READS=$OUTDIR/pooled.fastq

mkdir -p "$OUTDIR"
cd "$OUTDIR"        # all generated files land here

# 1) Targeted reference: one contig per gene (symbol-named), genomic interval w/ introns
awk '!seen[$4]++' "$BED" > panel.uniq.bed
: > kzfp_targets.fa
while read chr start end gene rest; do
    samtools faidx "$GENOME" "$chr:$((start+1))-$end" | sed "1s|.*|>$gene|" >> kzfp_targets.fa
done < panel.uniq.bed
samtools faidx kzfp_targets.fa

# 2) Splice-aware alignment, one placement per read (paralogs can't steal reads)
minimap2 -t "$THREADS" -ax splice -k14 --secondary=no -ub kzfp_targets.fa "$READS" \
  | samtools sort -@ "$THREADS" -o kzfp_isoforms.sorted.bam -
samtools index kzfp_isoforms.sorted.bam

# 3) Quick confirmation reads landed
echo "=== flagstat ==="
samtools flagstat kzfp_isoforms.sorted.bam
echo "=== reads per gene (top 20) ==="
samtools idxstats kzfp_isoforms.sorted.bam | sort -k3,3nr | head -20

echo "Done. In IGV: load kzfp_targets.fa as genome, then kzfp_isoforms.sorted.bam."