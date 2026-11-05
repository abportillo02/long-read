#!/bin/bash
#SBATCH --job-name=lift_grcm39_v2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/scripts/logs/lift_grcm39_v2_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
I="$BASE/genotype-reference/sample-refs/intermediate"
M="$BASE/mapping/methyl"
OUT="$M/grcm39_v2"                                   # CHANGED
G="$BASE/genotype-reference/GRCm39.fa"               # CHANGED (CG check)
LIBS="$BASE/mapping/libraries.tsv"

export PATH="$PATH:$HOME/miniconda3/bin"             # CHANGED (keeps active env first)

for t in liftOver chainSwap samtools; do             # CHANGED
    command -v "$t" >/dev/null || { echo "ERROR: $t not found on PATH" >&2; exit 1; }
done

mkdir -p "$OUT"

for f in \
    "$I/YRLP_0066.chain" "$I/YRLP_0067.chain" "$I/YRLP_0068.chain" \
    "$I/YRLP_0066.floxed_to_recomb.chain" "$G" "$G.fai" "$LIBS"; do
    [[ -s "$f" ]] || { echo "ERROR: Missing or empty file: $f" >&2; exit 1; }
done

# CHANGED: reverse the chains. The originals convert GRCm39 -> sample;
# the CpG coordinates are in sample space, so they need sample -> GRCm39.
C="$OUT/chains"; mkdir -p "$C"
for r in 0066 0067 0068; do
    chainSwap "$I/YRLP_${r}.chain" "$C/YRLP_${r}.to_grcm39.chain"
done
chainSwap "$I/YRLP_0066.floxed_to_recomb.chain" "$C/YRLP_0066.recomb_to_floxed.chain"

# CHANGED: percent of sampled sites reading CG in GRCm39
cg_pct() {
    local f=$1 reg=$2
    { if [[ $reg == chr15 ]]; then grep -w '^chr15' "$f" || true; else cat "$f"; fi; } \
        | shuf -n 50000 | awk '{print $1":"$2+1"-"$3}' > "$OUT/.cg_regions.$$"
    samtools faidx "$G" -r "$OUT/.cg_regions.$$" | grep -v '>' | tr a-z A-Z \
        | awk '{n++; if ($0=="CG") c++} END {if (n>0) printf "%.3f", 100*c/n; else printf "0"}'
    rm -f "$OUT/.cg_regions.$$"
}

SUMMARY="$OUT/liftover_summary.tsv"
printf 'library\tcpgs_in\tcpgs_uniquely_lifted\tpct_unique_lifted\tcg_pct_all\tcg_pct_chr15\n' > "$SUMMARY"   # CHANGED

while IFS=$'\t' read -r LIB ASSAY GROUP REFID; do
    [[ -n "${LIB:-}" ]] || continue
    echo "Processing $LIB ($ASSAY; $GROUP; $REFID)"

    INPUT="$M/${LIB}_CpG.bedGraph"
    [[ -s "$INPUT" ]] || { echo "ERROR: Missing input: $INPUT" >&2; exit 1; }

    awk 'BEGIN{OFS="\t"}
         $1 ~ /^chr([0-9]+|X|Y)$/ { print $1, $2, $3, NR "|" $5 ":" $6 }' \
         "$INPUT" > "$OUT/$LIB.sample.bed"

    N_IN=$(wc -l < "$OUT/$LIB.sample.bed")

    if [[ "$REFID" == "YRLP_0066" ]]; then
        awk '$1=="chr15"' "$OUT/$LIB.sample.bed" > "$OUT/$LIB.chr15.recomb.bed"
        awk '$1!="chr15"' "$OUT/$LIB.sample.bed" > "$OUT/$LIB.nonchr15.bed"

        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.nonchr15.bed" "$C/YRLP_0066.to_grcm39.chain" \
            "$OUT/$LIB.nonchr15.grcm39.multi.bed" "$OUT/$LIB.nonchr15.unmapped.bed"       # CHANGED chain

        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.chr15.recomb.bed" "$C/YRLP_0066.recomb_to_floxed.chain" \
            "$OUT/$LIB.chr15.floxed.multi.bed" "$OUT/$LIB.chr15.stage1.unmapped.bed"      # CHANGED chain

        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.chr15.floxed.multi.bed" "$C/YRLP_0066.to_grcm39.chain" \
            "$OUT/$LIB.chr15.grcm39.multi.bed" "$OUT/$LIB.chr15.stage2.unmapped.bed"      # CHANGED chain

        cat "$OUT/$LIB.nonchr15.grcm39.multi.bed" "$OUT/$LIB.chr15.grcm39.multi.bed" \
            > "$OUT/$LIB.all.grcm39.multi.bed"
    else
        liftOver -multiple -minMatch=1.0 \
            "$OUT/$LIB.sample.bed" "$C/${REFID}.to_grcm39.chain" \
            "$OUT/$LIB.all.grcm39.multi.bed" "$OUT/$LIB.unmapped.bed"                     # CHANGED chain
    fi

    awk 'BEGIN{FS=OFS="\t"} { split($4, a, /[|]/); n[a[1]]++ }
         END { for (id in n) if (n[id] == 1) print id }' \
         "$OUT/$LIB.all.grcm39.multi.bed" > "$OUT/$LIB.unique_ids.txt"

    awk 'BEGIN{FS=OFS="\t"}
         FILENAME == ARGV[1] { keep[$1] = 1; next }
         { split($4, a, /[|]/)
           if (keep[a[1]] && ($3-$2)==2) {
               split(a[2], c, ":")
               if (c[1] ~ /^[0-9]+$/ && c[2] ~ /^[0-9]+$/) print $1, $2, $3, c[1], c[2]
           } }' "$OUT/$LIB.unique_ids.txt" "$OUT/$LIB.all.grcm39.multi.bed" \
        | sort -k1,1 -k2,2n > "$OUT/$LIB.grcm39.counts.bed"

    N_OUT=$(wc -l < "$OUT/$LIB.grcm39.counts.bed")
    PCT=$(awk -v n="$N_OUT" -v d="$N_IN" 'BEGIN{ if (d>0) printf "%.4f", 100*n/d; else printf "NA" }')

    # CHANGED: coordinate check, stop the job if either region is not CG
    CG_ALL=$(cg_pct "$OUT/$LIB.grcm39.counts.bed" all)
    CG_15=$(cg_pct "$OUT/$LIB.grcm39.counts.bed" chr15)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$LIB" "$N_IN" "$N_OUT" "$PCT" "$CG_ALL" "$CG_15" >> "$SUMMARY"
    echo "  Lifted: $N_OUT of $N_IN ($PCT%); CG check all=$CG_ALL% chr15=$CG_15%"
    for v in "$CG_ALL" "$CG_15"; do
        awk -v x="$v" 'BEGIN{exit !(x>=99)}' || { echo "ERROR: CG check failed for $LIB" >&2; exit 1; }
    done

done < "$LIBS"

echo
echo "Finished. Summary: $SUMMARY"