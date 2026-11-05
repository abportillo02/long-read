#!/bin/bash
#SBATCH --job-name=bwameth_index
#SBATCH --array=0-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=BWAMETH_IDX_%A_%a.out
#SBATCH --error=BWAMETH_IDX_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

## After finalizing the reference FASTA for each sample, this script indexes it with bwa-meth and samtools, adding control sequences to the reference.

module load samtools
module load BWA/0.7.17-foss-2018b
module load Python
source ~/bwameth_env/bin/activate

set -euo pipefail

BWAMETH=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/20251029_LH00295_0299_B237232LT4/software/bwa-meth/bwameth.py
which python bwa samtools
python -c "import toolshed; print('toolshed ok')"
python $BWAMETH --help > /dev/null
module list 2>&1

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/sample-refs
CONTROLS=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/methyl_controls.fa     # EDIT: control FASTA(s), kit-supplied sequences

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

mkdir -p $BASE/mapping-refs
cd $BASE/mapping-refs

pip freeze > ${S}_bwameth_pip_freeze.txt

# Control contig names must not clash with any chr name
grep ">" $CONTROLS
if grep ">" $CONTROLS | cut -d' ' -f1 | tr -d '>' | grep -x -F -f <(cut -f1 $BASE/final/${S}_ref.fa.fai); then
  echo "FAILED: control contig name clashes with reference"; exit 1
fi

cat $BASE/final/${S}_ref.fa $CONTROLS > ${S}_ref.withcontrols.fa
samtools faidx ${S}_ref.withcontrols.fa
python $BWAMETH index ${S}_ref.withcontrols.fa

ls -lh ${S}_ref.withcontrols.fa*
echo "DONE"