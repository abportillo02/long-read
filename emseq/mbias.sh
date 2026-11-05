#!/bin/bash
#SBATCH --job-name=mbias
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=mbias_%j.out
#SBATCH --error=mbias_%j.err
#SBATCH	--mail-type=all
#SBATCH	--mail-user=abportillo@coh.org

set -euo pipefail
export PATH=/home/abportillo/.local/share/mamba/envs/methyldackel_env/bin:$PATH

EM=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REFDIR=$EM/genotype-reference/sample-refs/mapping-refs
TABLE=$EM/mapping/libraries.tsv
cd $EM/mapping
mkdir -p mbias

date
echo "MethylDackel: $(MethylDackel --version 2>&1 | head -n 1 || true)"

for L in $(cut -f1 $TABLE); do
  R=$(awk -F'\t' -v l="$L" '$1==l{print $4}' $TABLE)
  [ -s ${L}.markdup.bam ] || { echo "missing BAM for $L"; exit 1; }
  echo "mbias $L (reference $R)"
  MethylDackel mbias -r chr1 \
    $REFDIR/${R}_ref.withcontrols.fa ${L}.markdup.bam mbias/${L} \
    > mbias/${L}.mbias.txt 2> mbias/${L}.mbias.log
done
echo "DONE"