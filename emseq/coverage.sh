#!/bin/bash
#SBATCH --job-name=COV
#SBATCH --array=0-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/wgs-alignment/logs/COV_%A_%a.out
#SBATCH --error=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/wgs-alignment/logs/COV_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo 

sset -eo pipefail

source /home/abportillo/miniconda3/etc/profile.d/conda.sh
conda activate mamba_abner_BC

set -u

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
BAMDIR=$BASE/wgs-alignment
OUT=$BASE/coverage
mkdir -p $OUT logs

SAMPLES=(YRLP_0066 YRLP_0067 YRLP_0068)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Per-chromosome table (kept for reference)
samtools coverage $BAMDIR/${S}.markdup.bam > $OUT/${S}.coverage.txt

# Genome-wide mean depth, weighted by chromosome length, main chromosomes only
MEAN=$(awk 'NR>1 && $1 ~ /^(chr)?([0-9]+|X)$/ {sum+=$3*$7; len+=$3} END {print sum/len}' $OUT/${S}.coverage.txt)

echo -e "${S}\t${MEAN}" > $OUT/${S}.meandepth.txt
echo "${S} mean depth: ${MEAN}"

# check results 
# cat /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/coverage/YRLP_00*.meandepth.txt

# C