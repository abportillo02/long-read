#!/bin/bash
#SBATCH --job-name=offtarget_annotate
#SBATCH -p all
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=offtarget_annotate%j.out
#SBATCH --error=offtarget_annotate%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org


# annotate_offtarget.sh
# Rank the top read pile-ups in a genome BAM and label them with gene names.
# Tells you what got sequenced instead of (or alongside) your KZFP targets.
#
# Usage:
#   ./annotate_offtarget.sh <genome.sorted.bam> [min_depth] [n_top]
#
# Example:
#   ./annotate_offtarget.sh genome.sorted.bam 5 30
#
# Defaults: min_depth=5, n_top=20
#
# Run from inside the conda env that has bedtools + samtools:
#   conda activate mamba_abner_BC

set -euo pipefail

# ---- config (edit paths if they move) -------------------------------------
GTF="/home/abportillo/genomes/hg38/gencode.v43.annotation.gtf"
GENES_BED="/home/abportillo/github_repo/kzfp_long_read/results/kzfp_signal/gencode.v43.genes.bed"  # cache, built once
# ---------------------------------------------------------------------------

BAM="${1:?usage: annotate_offtarget.sh <genome.sorted.bam> [min_depth] [n_top]}"
MIN_DEPTH="${2:-5}"
N_TOP="${3:-20}"

[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }

# 1. Build the gene-level BED cache once (collapses GTF to ~60k lines).
if [[ ! -s "$GENES_BED" ]]; then
  echo ">> Building gene BED cache (one time): $GENES_BED" >&2
  [[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
  awk -F'\t' '$3=="gene"{
    name=$9; sub(/.*gene_name "/,"",name); sub(/".*/,"",name);
    print $1"\t"$4-1"\t"$5"\t"name
  }' "$GTF" > "$GENES_BED"
fi

# 2. chr-prefix sanity check (mismatch => empty results, not a real verdict).
bam_chr=$(samtools idxstats "$BAM" | head -1 | cut -f1)
bed_chr=$(head -1 "$GENES_BED" | cut -f1)
case "$bam_chr" in chr*) bam_pref=chr;; *) bam_pref=plain;; esac
case "$bed_chr" in chr*) bed_pref=chr;; *) bed_pref=plain;; esac
if [[ "$bam_pref" != "$bed_pref" ]]; then
  echo "WARNING: contig naming mismatch (BAM='$bam_chr' vs BED='$bed_chr')." >&2
  echo "         Results will be empty/wrong until both use the same convention." >&2
fi

# 3. Find covered intervals, keep the deepest, annotate with gene names.
echo ">> Top $N_TOP pile-ups (min depth $MIN_DEPTH) in $BAM" >&2
echo -e "region\tmax_depth\tgene"

bedtools genomecov -ibam "$BAM" -bga \
  | awk -v d="$MIN_DEPTH" '$4>=d' \
  | bedtools merge -d 50 -c 4 -o max \
  | sort -k4,4nr \
  | head -n "$N_TOP" \
  | bedtools intersect -a - -b "$GENES_BED" -wa -wb \
  | awk '{print $1":"$2"-"$3"\t"$4"\t"$8}' \
  | sort -t$'\t' -k2,2nr