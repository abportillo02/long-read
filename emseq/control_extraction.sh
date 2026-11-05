#!/bin/bash
#SBATCH --job-name=control_extraction
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=control_extraction_%j.out
#SBATCH --error=control_extraction_%j.err

set -euo pipefail
export PATH=/home/abportillo/.local/share/mamba/envs/methyldackel_env/bin:$PATH
which MethylDackel samtools

EM=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REFDIR=$EM/genotype-reference/sample-refs/mapping-refs
TABLE=$EM/mapping/libraries.tsv
cd $EM/mapping
mkdir -p controls

for L in $(cut -f1 $TABLE); do
  R=$(awk -F'\t' -v l="$L" '$1==l{print $4}' $TABLE)
  [ -n "$R" ] || { echo "no reference for $L"; exit 1; }
  [ -s ${L}.markdup.bam ] || { echo "missing BAM for $L"; exit 1; }
  for C in phage_lambda plasmid_puc19c phage_Xp12 phage_T4; do
    MethylDackel extract --minDepth 5 -r $C -o controls/${L}_${C} \
      $REFDIR/${R}_ref.withcontrols.fa ${L}.markdup.bam
  done
done

{
printf "library\ttable_label\tlambda\tpUC19\tXp12\tT4\n"
for L in $(cut -f1 $TABLE); do
  A=$(awk -F'\t' -v l="$L" '$1==l{print $2}' $TABLE)
  row="$L\t$A"
  for C in phage_lambda plasmid_puc19c phage_Xp12 phage_T4; do
    f=controls/${L}_${C}_CpG.bedGraph
    if [ -s $f ]; then
      v=$(awk 'NR>1{m+=$5; t+=$5+$6} END{if(t>0) printf "%.1f", 100*m/t; else print "NA"}' $f)
    else
      v=NA
    fi
    row="$row\t$v"
  done
  printf "$row\n"
done
} | tee controls/control_summary.tsv
echo "DONE"