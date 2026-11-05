#!/bin/bash
#SBATCH --job-name=bwameth_mapping
#SBATCH --array=0,6,7,8
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/MAP_%A_%a.out
#SBATCH --error=logs/MAP_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

module load samtools
module load BWA/0.7.17-foss-2018b
module load Python
source ~/bwameth_env/bin/activate

set -euo pipefail

BWAMETH=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/20251029_LH00295_0299_B237232LT4/software/bwa-meth/bwameth.py
EM=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
REFDIR=$EM/genotype-reference/sample-refs/mapping-refs
FQDIR=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/20251029_LH00295_0299_B237232LT4/raw-data/fastq
TABLE=$EM/mapping/libraries.tsv  
[ -s $TABLE ] || { echo "library table missing or empty: $TABLE"; exit 1; }        # columns: library_id, assay, strain, reference
OUT=$EM/mapping
mkdir -p $OUT

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TABLE)
L=$(echo "$LINE" | cut -f1)
R=$(echo "$LINE" | cut -f4)

R1=$(ls $FQDIR/${L}_*_R1_001.fastq.gz)
R2=$(ls $FQDIR/${L}_*_R2_001.fastq.gz)
[ "$(echo "$R1" | wc -l)" -eq 1 ] && [ "$(echo "$R2" | wc -l)" -eq 1 ] || { echo "expected one R1/R2 for $L"; exit 1; }
[ -n "$R" ] || { echo "no reference set for $L"; exit 1; }
[ -s $REFDIR/${R}_ref.withcontrols.fa.bwameth.c2t.sa ] || { echo "index missing for $R"; exit 1; }

echo "library $L ($(echo "$LINE" | cut -f2), $(echo "$LINE" | cut -f3)) -> reference $R"
which python bwa samtools
python -c "import toolshed; print('toolshed ok')"
module list 2>&1
pip freeze > $OUT/${L}_pip_freeze.txt

   python $BWAMETH --reference $REFDIR/${R}_ref.withcontrols.fa -t 12 \
--read-group "@RG\tID:${L}\tSM:${L}\tPL:ILLUMINA" $R1 $R2 \
| samtools fixmate -m -@ 2 - - \
| samtools sort -@ 3 -m 3G -T $OUT/${L}.tmp - \
| samtools markdup -@ 2 - $OUT/${L}.markdup.bam

samtools index $OUT/${L}.markdup.bam
samtools flagstat $OUT/${L}.markdup.bam > $OUT/${L}.flagstat.txt
samtools idxstats $OUT/${L}.markdup.bam > $OUT/${L}.idxstats.txt
echo "reads on control contigs:"
grep -E "^(phage|plasmid)" $OUT/${L}.idxstats.txt
echo "DONE"