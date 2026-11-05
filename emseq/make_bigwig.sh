#!/bin/bash
#SBATCH --job-name=BW
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org


# Make IGV-ready BigWig tracks from lifted (GRCm39) per-CpG methylation counts.
#
# For each library this writes two tracks:
#   <lib>_<assay>_<group>_methylation.bw   percent methylation (0-100) per CpG
#   <lib>_<assay>_<group>_depth.bw         read depth per CpG (reliability track)
# plus a .bigwig.params.txt file recording the inputs, parameters and tool version.
#
# Usage:
#   Slurm array (one task per line of mapping/libraries.tsv):
#       sbatch make_bigwigs.sh
#   Single library, interactive:
#       bash make_bigwigs.sh YRLP_0064
#
# Optional environment overrides: PROJECT, MIN_DEPTH (default 5), BG2BW (path to
# bedGraphToBigWig), COUNTS_DIR, COUNTS_SUFFIX.

set -euo pipefail

PROJECT="${PROJECT:-/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test}"
FASTA="$PROJECT/genotype-reference/GRCm39.fa"
LIBTSV="$PROJECT/mapping/libraries.tsv"          # no header: library, assay, group, reference
COUNTS_DIR="${COUNTS_DIR:-$PROJECT/mapping/methyl/grcm39_v2}"
COUNTS_SUFFIX="${COUNTS_SUFFIX:-.grcm39.counts.bed}"   # columns: chrom start end meth unmeth
OUTDIR="$COUNTS_DIR/bigwig"
SIZES="$OUTDIR/GRCm39.chrom.sizes"
MIN_DEPTH="${MIN_DEPTH:-5}"
BG2BW="${BG2BW:-bedGraphToBigWig}"
MAX_BAD_FRACTION="${MAX_BAD_FRACTION:-0.001}"   # abort if more than this fraction of CpGs are dropped

# ---- which library --------------------------------------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    LIB="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIBTSV" | cut -f1)"
else
    LIB="${1:?usage: bash make_bigwigs.sh YRLP_00NN}"
fi
[[ "$LIB" == YRLP_* ]] || LIB="YRLP_$LIB"

read -r ASSAY GROUP < <(
    awk -F'\t' -v l="$LIB" '($1 == l || "YRLP_" $1 == l) { print $2, $3; exit }' "$LIBTSV"
)
[[ -n "${ASSAY:-}" ]] || { echo "ERROR: $LIB not found in $LIBTSV" >&2; exit 1; }
SAFE_NAME="$(printf '%s_%s_%s' "$LIB" "$ASSAY" "$GROUP" | tr -c 'A-Za-z0-9._\n-' '_')"

# ---- checks ---------------------------------------------------------------
IN="$COUNTS_DIR/${LIB}${COUNTS_SUFFIX}"
[[ -s "$IN" ]]    || { echo "ERROR: counts file not found: $IN" >&2; exit 1; }
[[ -s "$SIZES" ]] || {
    echo "ERROR: $SIZES missing. Create it once with:" >&2
    echo "  mkdir -p $OUTDIR; module load samtools" >&2
    echo "  [ -s $FASTA.fai ] || samtools faidx $FASTA" >&2
    echo "  cut -f1,2 $FASTA.fai > $SIZES" >&2
    exit 1
}
command -v "$BG2BW" >/dev/null || { echo "ERROR: $BG2BW not found on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/bw_${LIB}.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

echo "Library $LIB ($ASSAY, $GROUP); min depth $MIN_DEPTH; input $IN"

# ---- step 1: depth filter, percent methylation, drop unknown chromosomes ---
# bigWig stores one value per interval, so each CpG becomes a 1-bp interval with
# percent = 100 * meth / (meth + unmeth). Chromosomes absent from the FASTA index
# (or intervals past the chromosome end) cannot be written and are counted.
awk -F'\t' -v OFS='\t' -v min="$MIN_DEPTH" -v sizes="$SIZES" -v stats="$TMP/stats1" '
BEGIN {
    while ((getline line < sizes) > 0) { split(line, a, "\t"); known[a[1]] = a[2] + 0 }
}
$2 !~ /^[0-9]+$/ { next }                      # skip any header or track line
{
    n++
    if (!($1 in known) || ($3 + 0) > known[$1]) { bad++; next }
    d = $4 + $5
    if (d >= min) { kept++; printf "%s\t%d\t%d\t%.4f\t%d\n", $1, $2, $3, 100 * $4 / d, d }
}
END { print n + 0, bad + 0, kept + 0 > stats }
' "$IN" | LC_ALL=C sort -S 4G -T "$TMP" -k1,1 -k2,2n > "$TMP/sorted.tsv"

read -r N_TOTAL N_BAD N_KEPT < "$TMP/stats1"
awk -v n="$N_TOTAL" -v b="$N_BAD" -v max="$MAX_BAD_FRACTION" \
    'BEGIN { exit !(n > 0 && b / n <= max) }' || {
    echo "ERROR: $N_BAD of $N_TOTAL CpGs dropped (unknown chromosome or past chromosome end)." >&2
    echo "Check that chromosome names in $IN match $FASTA (chr1 vs 1)." >&2
    exit 1
}

# ---- step 2: remove overlapping or duplicate intervals ---------------------
# bedGraphToBigWig refuses overlapping intervals; keep the first of any overlap.
awk -F'\t' -v OFS='\t' -v stats="$TMP/stats2" '
{ if ($1 == pc && $2 < pe) { dup++; next } pc = $1; pe = $3; print }
END { print dup + 0 > stats }
' "$TMP/sorted.tsv" > "$TMP/final.tsv"
read -r N_DUP < "$TMP/stats2"

# ---- step 3: write bedGraphs and convert -----------------------------------
cut -f1-4 "$TMP/final.tsv"  > "$TMP/methylation.bedGraph"
cut -f1-3,5 "$TMP/final.tsv" > "$TMP/depth.bedGraph"

for track in methylation depth; do
    "$BG2BW" "$TMP/$track.bedGraph" "$SIZES" "$OUTDIR/${SAFE_NAME}_${track}.bw.tmp"
    mv "$OUTDIR/${SAFE_NAME}_${track}.bw.tmp" "$OUTDIR/${SAFE_NAME}_${track}.bw"
done

# ---- step 4: record how the tracks were made -------------------------------
BW_VERSION="$("$BG2BW" 2>&1 | head -1 || true)"
{
    echo "date:            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host:            $(hostname)"
    echo "script:          $(readlink -f "$0")"
    echo "library:         $LIB"
    echo "assay/group:     $ASSAY / $GROUP"
    echo "input:           $IN"
    echo "input_size_mtime: $(stat -c '%s bytes, %y' "$IN")"
    echo "chrom_sizes:     $SIZES (from $FASTA.fai)"
    echo "min_depth:       $MIN_DEPTH"
    echo "value:           100 * methylated / (methylated + unmethylated), per CpG"
    echo "tool:            $BW_VERSION"
    echo "cpgs_in_input:   $N_TOTAL"
    echo "cpgs_dropped_unknown_chrom_or_off_end: $N_BAD"
    echo "cpgs_passing_depth: $N_KEPT"
    echo "cpgs_dropped_overlap: $N_DUP"
    echo "outputs:         ${SAFE_NAME}_methylation.bw, ${SAFE_NAME}_depth.bw"
} > "$OUTDIR/${SAFE_NAME}.bigwig.params.txt"

echo "Done: $OUTDIR/${SAFE_NAME}_{methylation,depth}.bw ($N_KEPT CpGs; $N_DUP overlaps dropped)"