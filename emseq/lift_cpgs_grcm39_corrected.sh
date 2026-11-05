#!/bin/bash
#SBATCH --job-name=lift_grcm39_fix
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/scripts/logs/lift_grcm39_corrected_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
I="$BASE/genotype-reference/sample-refs/intermediate"
M="$BASE/mapping/methyl"
OUT="$M/grcm39_corrected"
LIBS="$BASE/mapping/libraries.tsv"

export PATH="$HOME/miniconda3/bin:$PATH"

command -v liftOver >/dev/null || {
    echo "ERROR: liftOver not found on PATH" >&2
    exit 1
}

mkdir -p "$OUT"

# Validate inputs before starting.
for f in \
    "$I/YRLP_0066.chain" \
    "$I/YRLP_0067.chain" \
    "$I/YRLP_0068.chain" \
    "$I/YRLP_0066.floxed_to_recomb.chain" \
    "$LIBS"; do
    [[ -s "$f" ]] || {
        echo "ERROR: Missing or empty file: $f" >&2
        exit 1
    }
done

SUMMARY="$OUT/liftover_summary.tsv"
printf 'library\tcpgs_in\tcpgs_uniquely_lifted\tpct_unique_lifted\n' \
    > "$SUMMARY"

while IFS=$'\t' read -r LIB ASSAY GROUP REFID; do
    [[ -n "${LIB:-}" ]] || continue

    echo "Processing $LIB ($ASSAY; $GROUP; $REFID)"

    INPUT="$M/${LIB}_CpG.bedGraph"
    [[ -s "$INPUT" ]] || {
        echo "ERROR: Missing input: $INPUT" >&2
        exit 1
    }

    # Retain standard chromosomes and carry a unique source ID plus counts.
    # This assumes methylated and unmethylated counts are columns 5 and 6,
    # as in the original script.
    awk 'BEGIN{OFS="\t"}
         $1 ~ /^chr([0-9]+|X|Y)$/ {
             print $1, $2, $3, NR "|" $5 ":" $6
         }' "$INPUT" > "$OUT/$LIB.sample.bed"

    N_IN=$(wc -l < "$OUT/$LIB.sample.bed")

    if [[ "$REFID" == "YRLP_0066" ]]; then
        # The recombined reference differs from the floxed reference
        # only on chr15.
        awk '$1=="chr15"' "$OUT/$LIB.sample.bed" \
            > "$OUT/$LIB.chr15.recomb.bed"

        awk '$1!="chr15"' "$OUT/$LIB.sample.bed" \
            > "$OUT/$LIB.nonchr15.bed"

        # Non-chr15 coordinates are unchanged by the recombination.
        # Original chain: sample/floxed query -> GRCm39 target.
        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.nonchr15.bed" \
            "$I/YRLP_0066.chain" \
            "$OUT/$LIB.nonchr15.grcm39.multi.bed" \
            "$OUT/$LIB.nonchr15.unmapped.bed"

        # Chr15 stage 1: recombined -> floxed.
        # The original special chain already has this direction.
        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.chr15.recomb.bed" \
            "$I/YRLP_0066.floxed_to_recomb.chain" \
            "$OUT/$LIB.chr15.floxed.multi.bed" \
            "$OUT/$LIB.chr15.stage1.unmapped.bed"

        # Chr15 stage 2: floxed -> GRCm39.
        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.chr15.floxed.multi.bed" \
            "$I/YRLP_0066.chain" \
            "$OUT/$LIB.chr15.grcm39.multi.bed" \
            "$OUT/$LIB.chr15.stage2.unmapped.bed"

        cat "$OUT/$LIB.nonchr15.grcm39.multi.bed" \
            "$OUT/$LIB.chr15.grcm39.multi.bed" \
            > "$OUT/$LIB.all.grcm39.multi.bed"

    else
        # Original chain maps the sample reference to GRCm39.
        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.sample.bed" \
            "$I/${REFID}.chain" \
            "$OUT/$LIB.all.grcm39.multi.bed" \
            "$OUT/$LIB.unmapped.bed"
    fi

    # Retain source IDs with exactly one final mapping.
    awk 'BEGIN{FS=OFS="\t"}
         {
             split($4, a, /[|]/)
             n[a[1]]++
         }
         END {
             for (id in n)
                 if (n[id] == 1) print id
         }' "$OUT/$LIB.all.grcm39.multi.bed" \
        > "$OUT/$LIB.unique_ids.txt"

    # Restore the two count columns, require a 2-bp interval,
    # and sort in GRCm39 coordinate order.
    awk 'BEGIN{FS=OFS="\t"}
         FILENAME == ARGV[1] {
             keep[$1] = 1
             next
         }
         {
             split($4, a, /[|]/)
             if (keep[a[1]] && ($3-$2)==2) {
                 split(a[2], c, ":")
                 if (c[1] ~ /^[0-9]+$/ && c[2] ~ /^[0-9]+$/)
                     print $1, $2, $3, c[1], c[2]
             }
         }' "$OUT/$LIB.unique_ids.txt" \
            "$OUT/$LIB.all.grcm39.multi.bed" \
        | sort -k1,1 -k2,2n \
        > "$OUT/$LIB.grcm39.counts.bed"

    N_OUT=$(wc -l < "$OUT/$LIB.grcm39.counts.bed")

    if (( N_IN > 0 )); then
        PCT=$(awk -v n="$N_OUT" -v d="$N_IN" \
            'BEGIN{printf "%.4f", 100*n/d}')
    else
        PCT="NA"
    fi

    printf '%s\t%s\t%s\t%s\n' \
        "$LIB" "$N_IN" "$N_OUT" "$PCT" >> "$SUMMARY"

    echo "  Input CpGs: $N_IN"
    echo "  Uniquely lifted CpGs: $N_OUT ($PCT%)"

done < "$LIBS"

echo
echo "Finished. Summary: $SUMMARY"
echo "Corrected bedGraphs: $OUT/*.grcm39.counts.bed"

