#!/bin/bash
#SBATCH --job-name=lift_grcm39
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/lift_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH="$HOME/miniconda3/bin:$PATH"
command -v liftOver chainSwap >/dev/null || { echo "liftOver/chainSwap not on PATH" >&2; exit 1; }

I=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/sample-refs/intermediate
MAP=../mapping; M=$MAP/methyl; OUT=$M/grcm39; C=$OUT/chains
mkdir -p "$C"

# Reverse the chains: GRCm39 -> sample becomes sample -> GRCm39
for r in 0066 0067 0068; do
  chainSwap $I/YRLP_$r.chain $C/YRLP_$r.sample_to_grcm39.chain
done
# 0066 only: final (Cre-recombined) -> floxed intermediate
chainSwap $I/YRLP_0066.floxed_to_recomb.chain $C/YRLP_0066.recomb_to_floxed.chain

{
  echo "date: $(date -Is)"
  echo "chains: $I"
  echo "liftOver options: -minMatch=1.0 (both bases of each CpG must map, one-to-one)"
  echo "tools:"; conda list 2>/dev/null | grep -E "ucsc-(liftover|chainswap)" || echo "  (conda list unavailable)"
} > $OUT/lift_params.txt
printf "library\tcpgs_in\tcpgs_lifted\tpct_lifted\n" > $OUT/lift_summary.tsv

while IFS=$'\t' read -r LIB ASSAY GROUP REFID <&3; do
  # standard chromosomes only (control contigs have no chain); counts ride in the name field
  awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print $1,$2,$3,$5":"$6}' \
    $M/${LIB}_CpG.bedGraph > $OUT/$LIB.sample.bed

  if [[ $REFID == YRLP_0066 ]]; then
    liftOver -minMatch=1.0 $OUT/$LIB.sample.bed $C/YRLP_0066.recomb_to_floxed.chain \
      $OUT/$LIB.floxed.bed $OUT/$LIB.unmapped1.bed
    liftOver -minMatch=1.0 $OUT/$LIB.floxed.bed $C/YRLP_0066.sample_to_grcm39.chain \
      $OUT/$LIB.lifted.bed $OUT/$LIB.unmapped2.bed
  else
    liftOver -minMatch=1.0 $OUT/$LIB.sample.bed $C/${REFID}.sample_to_grcm39.chain \
      $OUT/$LIB.lifted.bed $OUT/$LIB.unmapped2.bed
  fi

  # keep 2 bp CpGs only; columns: chr, start, end, methylated, unmethylated
  awk 'BEGIN{OFS="\t"} ($3-$2)==2 {split($4,a,":"); print $1,$2,$3,a[1],a[2]}' $OUT/$LIB.lifted.bed \
    | sort -k1,1 -k2,2n > $OUT/$LIB.grcm39.counts.bed

  n_in=$(wc -l < $OUT/$LIB.sample.bed); n_out=$(wc -l < $OUT/$LIB.grcm39.counts.bed)
  printf "%s\t%d\t%d\t%s\n" $LIB $n_in $n_out \
    "$(awk -v a=$n_in -v b=$n_out 'BEGIN{printf "%.2f", 100*b/a}')" >> $OUT/lift_summary.tsv
done 3< $MAP/libraries.tsv

cat $OUT/lift_summary.tsv