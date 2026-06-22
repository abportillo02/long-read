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
# Input : calls.bam from the gemini Dorado step (this run: LSK114 finishing on
#         PCB114-barcoded cDNA, capture-enriched; min-qscore 9).
#
# Goal 1 (did capture work?): per-gene read counts + on-target %.
#         -> TRUST THE GENOME ON-TARGET FRACTION, not the raw panel number.
#            Panel-FASTA mapping can be inflated by repeat/paralog mis-mapping
#            of KZFP zinc-finger arrays at low MAPQ; the -q 1 filter and the
#            genome cross-check guard against calling that "signal".
# Goal 2 (isoforms?): genome splice BAM for IGV, then pychopper+IsoQuant
#         (see NEXT STEP block). Preset must be PCS114 for SQK-PCB114.24.
#
# HOW TO READ SUMMARY.txt:
#   1. Genome on-target fraction  -> honest capture-efficiency number
#   2. genes_with_reads / ge20    -> breadth: how many of 117 genes got covered
#   3. Panel on-target (-q 1)      -> if it barely drops vs no filter, signal is real;
#                                     if it collapses, much of it was repeat mis-mapping
#
# Needs samtools + minimap2 + gawk (present in mamba_abner_BC).
# =============================================================================

set -uo pipefail

# --- Environment ---
set +u                                              # /etc/bashrc reads PS1, unset in batch jobs
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC
set -u

THREADS=${SLURM_CPUS_PER_TASK:-16}

# --- Paths ---
# NOTE: OG_Path (this run, holds calls.bam) is intentionally separate from the
# shared reference docs under long_read/. Confirm BOTH resolve before submit:
#   ls -d /home/abportillo/github_repo/*long*read*/
#   ls "$GENOME" "$PANEL_FASTA" "$GENE_LIST" "$GTF" "$CALLS_BAM"
OG_Path=/home/abportillo/github_repo/kzfp_long_read
GENOME=/home/abportillo/github_repo/long_read/docs/hg38_p14.fa
GTF=/home/abportillo/github_repo/long_read/docs/gencode.v43.annotation.gtf
PANEL_FASTA=/home/abportillo/github_repo/long_read/docs/primate_kzfp_ordered.fasta   # master-transcript panel, 117 genes, exon-concatenated
GENE_LIST=/home/abportillo/github_repo/long_read/docs/priority_kzfps_ordered.txt
OUTDIR=$OG_Path/results/kzfp_signal
CALLS_BAM=$OG_Path/calls.bam

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# --- Fail loudly on a bad path BEFORE minimap2 produces a silent broken BAM ---
if [[ ! -s "$CALLS_BAM" ]]; then
    echo "ERROR: $CALLS_BAM not found. Copy it from gemini or fix CALLS_BAM." >&2
    exit 1
fi
for f in "$GENOME" "$PANEL_FASTA" "$GENE_LIST" "$GTF"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty reference: $f" >&2; exit 1; }
done
# Sanity check the input BAM isn't truncated from transfer
samtools quickcheck "$CALLS_BAM" || { echo "ERROR: $CALLS_BAM failed quickcheck (truncated?)." >&2; exit 1; }

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

# Per-gene counts: PRIMARY reads only (-F 0x904) and MAPQ>=1 (-q 1) so that
# supplementary/repeat mis-mapped paralog reads are not counted as signal.
# Gene lengths joined from idxstats; genes with 0 reads are retained.
{
  echo -e "gene\tref_len\tprimary_reads"
  awk 'BEGIN{OFS="\t"}
       FNR==NR { if($1!="*"){len[$1]=$2; cnt[$1]=0} ; next }
       { cnt[$1]++ }
       END { for(g in len) print g, len[g], cnt[g] }' \
       <(samtools idxstats "$OUTDIR/panel.sorted.bam") \
       <(samtools view -F 0x904 -q 1 "$OUTDIR/panel.sorted.bam" | cut -f3) \
    | sort -t$'\t' -k3,3nr
} > "$OUTDIR/per_gene_read_counts.tsv"

panel_mapped=$(samtools view -c -F 0x904 -q 1 "$OUTDIR/panel.sorted.bam")
panel_mapped_nofilter=$(samtools view -c -F 0x904 "$OUTDIR/panel.sorted.bam")   # for the collapse check
genes_with_reads=$(awk -F'\t' 'NR>1 && $3>0'  "$OUTDIR/per_gene_read_counts.tsv" | wc -l)
genes_ge20=$(awk -F'\t' 'NR>1 && $3>=20' "$OUTDIR/per_gene_read_counts.tsv" | wc -l)

# ---------------------------------------------------------------------------
# 3. ALIGNMENT B -- genome, splice-aware (Goal 2: isoforms in IGV)
#    -ub (both strands); reads here are NOT pychopper-oriented, so both-strand
#    is correct. Optional junction guidance if paftools.js + k8 are available.
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
# 4. Genomic panel BED + genome on-target/coverage (the honest capture metric)
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
  echo "Total called reads:     $total_reads   (min-qscore 9; not comparable to q3 pilot yield)"
  echo
  echo "--- Panel master-transcript alignment (Goal 1) ---"
  echo "Reads mapped to panel (-q 1):     $panel_mapped"
  echo "Reads mapped to panel (no MAPQ):  $panel_mapped_nofilter   (if >> -q 1 value, signal is repeat-inflated)"
  [[ "$total_reads" -gt 0 ]] && \
    awk -v a="$panel_mapped" -v b="$total_reads" 'BEGIN{printf "On-target fraction (-q 1):        %.4f\n", a/b}'
  echo "Panel genes with >=1 read:        $genes_with_reads"
  echo "Panel genes with >=20 reads:      $genes_ge20"
  echo
  echo "--- Genome splice alignment (HONEST capture metric) ---"
  echo "Reads mapped to genome:           $genome_mapped"
  echo "Reads on panel (genome):          $genome_ontarget"
  [[ "$genome_mapped" -gt 0 ]] && \
    awk -v a="$genome_ontarget" -v b="$genome_mapped" 'BEGIN{printf "On-target (genome):               %.4f   <- capture verdict\n", a/b}'
  echo
  echo "Top genes by primary read count:"
  head -n 16 "$OUTDIR/per_gene_read_counts.tsv"
  echo "==========================================================="
} | tee "$OUTDIR/SUMMARY.txt"

echo "[$(date)] Done."
echo "Key outputs:"
echo "  $OUTDIR/SUMMARY.txt                 <- read this first (capture verdict)"
echo "  $OUTDIR/per_gene_read_counts.tsv    <- per-gene signal (primary, -q 1)"
echo "  $OUTDIR/panel.genes.coverage.tsv    <- per-gene genome coverage"
echo "  $OUTDIR/genome.sorted.bam (+ .bai)  <- load in IGV with hg38_p14 + GTF for isoforms"
echo
echo "  IGV note: load genome.sorted.bam against hg38_p14 (chr-prefixed names must match)."
echo "  If loading the FASTA itself as reference, run: samtools faidx $GENOME"

# ---------------------------------------------------------------------------
# NEXT STEP (only if signal is good) -- isoform calling on well-covered genes.
#
#   Library = PCB114 cDNA (SQK-PCB114.24), capture-enriched, LSK114-finished.
#   The PCB114 cDNA primers are on the INSIDE of the molecule, so pychopper
#   applies. PRESET MUST BE PCS114 (matches PCB114). A wrong preset finds ~no
#   full-length reads and produces garbage isoforms.
#
#   # 1. Orient/trim full-length cDNA. Read stats.txt: a healthy full-length %
#   #    confirms the PCB primers survived barcode trimming. If full-length is
#   #    near zero, the primers were trimmed off -> SKIP pychopper and run
#   #    IsoQuant on the unstranded genome.sorted.bam instead (lose orientation).
#   pychopper -k PCS114 -r "$OUTDIR/pychopper_report.pdf" -S "$OUTDIR/pychopper_stats.txt" \
#       "$OUTDIR/pooled.fastq" "$OUTDIR/full_length.fastq"
#   cat "$OUTDIR/pychopper_stats.txt"
#
#   # 2. Re-align oriented reads (-uf: forward-stranded after pychopper).
#   minimap2 -t "$THREADS" -ax splice -k14 -uf ${JUNCBED:+--junc-bed "$JUNCBED"} \
#       "$GENOME" "$OUTDIR/full_length.fastq" \
#     | samtools sort -@ "$THREADS" -o "$OUTDIR/fl.sorted.bam" -
#   samtools index "$OUTDIR/fl.sorted.bam"
#
#   # 3. Isoform calling. Only meaningful on genes that cleared a real coverage
#   #    bar (e.g. the genes_ge20 set); isoforms from 3-read genes are noise.
#   isoquant.py --reference "$GENOME" --genedb "$GTF" --data_type nanopore \
#       --bam "$OUTDIR/fl.sorted.bam" -o "$OUTDIR/isoquant"
# ---------------------------------------------------------------------------