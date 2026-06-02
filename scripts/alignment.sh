#!/bin/bash
#SBATCH --job-name=kzfp_align
#SBATCH -p all         
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=kzfp_align_%j.out
#SBATCH --error=kzfp_align_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

# =============================================================================
# STEP 2 of 2 -- ALIGNMENT + PILOT QC (CPU, run on apollo)
#
# Input : calls.bam from the gemini basecalling step.
# Output: per-gene read counts + on-target % (Goal 1: is there signal?)
#         genome splice BAM for IGV          (Goal 2: isoforms)
#
# Needs samtools + minimap2 + gawk (present in mamba_abner_BC).
# =============================================================================

set -uo pipefail

# --- Environment ---
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC   # has samtools + minimap2

THREADS=${SLURM_CPUS_PER_TASK:-16}

# --- Paths (EDIT to match apollo) ---
BASE=/home/abportillo/github_repo/long-read
GENOME=$BASE/docs/hg38_p14.fa
GTF=$BASE/docs/gencode.v43.annotation.gtf
PANEL_FASTA=$BASE/docs/primate_kzfp_probes.fasta     #  master-transcript panel
GENE_LIST=$BASE/docs/priority_kzfps_ordered.txt
OUTDIR=$BASE/results/pilot_kzfp
CALLS_BAM=$OUTDIR/calls.bam                         

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ ! -s "$CALLS_BAM" ]]; then
    echo "ERROR: $CALLS_BAM not found. Copy it from gemini or fix CALLS_BAM." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 1. BAM -> pooled fastq
# ---------------------------------------------------------------------------
echo "[$(date)] Extracting reads from $CALLS_BAM..."
samtools fastq -@ "$THREADS" "$CALLS_BAM" > "$OUTDIR/pooled.fastq"
total_reads=$(( $(wc -l < "$OUTDIR/pooled.fastq") / 4 ))
echo "[$(date)] Pooled reads: $total_reads"
[[ "$total_reads" -eq 0 ]] && { echo "ERROR: no reads in $CALLS_BAM" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 2. ALIGNMENT A -- panel master-transcript reference (Goal 1: signal)
#    map-ont (panel seqs are exon-concatenated; no genomic introns).
#    --secondary=no so each read is counted once, by best-matching gene.
# ---------------------------------------------------------------------------
echo "[$(date)] Aligning to panel master-transcripts..."
minimap2 -t "$THREADS" -ax map-ont --secondary=no "$PANEL_FASTA" "$OUTDIR/pooled.fastq" \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/panel.sorted.bam" -
samtools index "$OUTDIR/panel.sorted.bam"

samtools idxstats "$OUTDIR/panel.sorted.bam" \
  | awk 'BEGIN{OFS="\t"; print "gene","ref_len","mapped_reads","unmapped"} $1!="*"{print $1,$2,$3,$4}' \
  | sort -t$'\t' -k3,3nr > "$OUTDIR/per_gene_read_counts.tsv"

panel_mapped=$(samtools view -c -F 0x904 "$OUTDIR/panel.sorted.bam")
genes_with_reads=$(awk -F'\t' 'NR>1 && $3>0'  "$OUTDIR/per_gene_read_counts.tsv" | wc -l)
genes_ge20=$(awk -F'\t' 'NR>1 && $3>=20' "$OUTDIR/per_gene_read_counts.tsv" | wc -l)

# ---------------------------------------------------------------------------
# 3. ALIGNMENT B -- genome, splice-aware (Goal 2: isoforms in IGV)
#    -ub (both strands); reads are not orientation-corrected here.
#    Optional junction guidance if paftools.js + k8 are available.
# ---------------------------------------------------------------------------
JUNCBED="$OUTDIR/gencode.junctions.bed"
if command -v paftools.js >/dev/null 2>&1 && command -v k8 >/dev/null 2>&1; then
    echo "[$(date)] Building junction BED from GTF..."
    paftools.js gff2bed "$GTF" > "$JUNCBED" 2>/dev/null || JUNCBED=""
else
    JUNCBED=""
fi

echo "[$(date)] Aligning to genome (splice)..."
minimap2 -t "$THREADS" -ax splice -k14 -ub ${JUNCBED:+--junc-bed "$JUNCBED"} \
    "$GENOME" "$OUTDIR/pooled.fastq" \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/genome.sorted.bam" -
samtools index "$OUTDIR/genome.sorted.bam"

# ---------------------------------------------------------------------------
# 4. Genomic panel BED + genome-based on-target/coverage (cross-check)
# ---------------------------------------------------------------------------
echo "[$(date)] Building genomic panel BED..."
tr -d '\r' < "$GENE_LIST" \
  | awk '{gsub(/^[ \t]+|[ \t]+$/,"",$0)} $0!=""{print "gene_name \"" $0 "\""}' \
  > "$OUTDIR/patterns.txt"

grep -Ff "$OUTDIR/patterns.txt" "$GTF" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{
        if (match($9,/gene_name "([^"]+)"/,g) && g[1]!="")
            print $1, $4-1, $5, g[1], ".", $7
    }' \
  | sort -k1,1 -k2,2n > "$OUTDIR/panel.genes.bed"

genome_mapped=$(samtools view -c -F 0x904 "$OUTDIR/genome.sorted.bam")
genome_ontarget=0
if [[ -s "$OUTDIR/panel.genes.bed" ]]; then
    genome_ontarget=$(samtools view -c -F 0x904 -L "$OUTDIR/panel.genes.bed" "$OUTDIR/genome.sorted.bam")
    samtools bedcov "$OUTDIR/panel.genes.bed" "$OUTDIR/genome.sorted.bam" \
      | awk 'BEGIN{OFS="\t"; print "chrom","start","end","gene","score","strand","summed_depth"}{print}' \
      > "$OUTDIR/panel.genes.coverage.tsv"
fi

# ---------------------------------------------------------------------------
# 5. Summary
# ---------------------------------------------------------------------------
{
  echo "==================== KZFP PILOT SUMMARY ===================="
  echo "Date:                   $(date)"
  echo "Total basecalled reads: $total_reads"
  echo
  echo "--- Panel master-transcript alignment (Goal 1) ---"
  echo "Reads mapped to panel:  $panel_mapped"
  [[ "$total_reads" -gt 0 ]] && \
    awk -v a="$panel_mapped" -v b="$total_reads" 'BEGIN{printf "On-target fraction:     %.4f\n", a/b}'
  echo "Panel genes with >=1 read:   $genes_with_reads"
  echo "Panel genes with >=20 reads: $genes_ge20"
  echo
  echo "--- Genome splice alignment (cross-check) ---"
  echo "Reads mapped to genome: $genome_mapped"
  echo "Reads on panel (genome):$genome_ontarget"
  [[ "$genome_mapped" -gt 0 ]] && \
    awk -v a="$genome_ontarget" -v b="$genome_mapped" 'BEGIN{printf "On-target (genome):     %.4f\n", a/b}'
  echo
  echo "Top genes by read count:"
  head -n 16 "$OUTDIR/per_gene_read_counts.tsv"
  echo "==========================================================="
} | tee "$OUTDIR/SUMMARY.txt"

echo "[$(date)] Done."
echo "Key outputs:"
echo "  $OUTDIR/SUMMARY.txt                 <- read this first"
echo "  $OUTDIR/per_gene_read_counts.tsv    <- per-gene signal"
echo "  $OUTDIR/panel.genes.coverage.tsv    <- per-gene genome coverage"
echo "  $OUTDIR/genome.sorted.bam (+ .bai)  <- load in IGV with the GTF for isoforms"

# ---------------------------------------------------------------------------
# NEXT STEP (only if signal is good) -- isoform calling on covered genes.
#   pychopper -k PCS111 pooled.fastq full_length.fastq      # verify primer preset for PCB/PCS114
#   minimap2 -ax splice -k14 -uf ${JUNCBED:+--junc-bed "$JUNCBED"} \
#       "$GENOME" full_length.fastq | samtools sort -o fl.sorted.bam - ; samtools index fl.sorted.bam
#   isoquant.py --reference "$GENOME" --genedb "$GTF" --data_type nanopore \
#       --bam fl.sorted.bam -o "$OUTDIR/isoquant"
# ---------------------------------------------------------------------------