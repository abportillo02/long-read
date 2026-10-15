#!/bin/bash
#SBATCH --job-name=control_qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=control_qc_%j.out
#SBATCH --error=control_qc_%j.err

set -euo pipefail
export PATH=/home/abportillo/.local/share/mamba/envs/methyldackel_env/bin:$PATH

EM=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REFDIR=$EM/genotype-reference/sample-refs/mapping-refs
TABLE=$EM/mapping/libraries.tsv
cd $EM/mapping
mkdir -p controls

date
echo "MethylDackel: $(MethylDackel --version 2>&1 | head -n 1 || true)"
samtools --version | head -n 1

OUTTSV=controls/control_qc.tsv
printf "library\tassay\tcontrol\tcontext\tsites\tmethylated\tunmethylated\ttotal\tpct_methylated\n" > $OUTTSV

for L in $(cut -f1 $TABLE); do
  R=$(awk -F'\t' -v l="$L" '$1==l{print $4}' $TABLE)
  A=$(awk -F'\t' -v l="$L" '$1==l{print $2}' $TABLE)
  [ -n "$R" ] || { echo "no reference for $L"; exit 1; }
  [ -s ${L}.markdup.bam ] || { echo "missing BAM for $L"; exit 1; }
  for C in phage_lambda plasmid_puc19c phage_Xp12 phage_T4; do
    MethylDackel extract --CHG --CHH --minDepth 1 -r $C -o controls/${L}_${C} \
      $REFDIR/${R}_ref.withcontrols.fa ${L}.markdup.bam
    for X in CpG CHG CHH; do
      f=controls/${L}_${C}_${X}.bedGraph
      if [ -s "$f" ]; then
        awk -v l="$L" -v a="$A" -v c="$C" -v x="$X" \
          'NR>1{m+=$5; u+=$6; s++} END{t=m+u; printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n", l,a,c,x,s,m,u,t,(t>0 ? sprintf("%.2f",100*m/t) : "NA")}' "$f" >> $OUTTSV
      else
        printf "%s\t%s\t%s\t%s\t0\t0\t0\t0\tNA\n" "$L" "$A" "$C" "$X" >> $OUTTSV
      fi
    done
  done
done
echo "DONE"