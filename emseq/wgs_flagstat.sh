#!/bin/bash
#SBATCH --job-name=WGS_flagstat
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=wgs_flagstat_%j.out
#SBATCH --error=wgs_flagstat_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mamba_abner_BC

BAM_DIR="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/wgs-alignment"

cd "$BAM_DIR"

for sample in YRLP_0066 YRLP_0067 YRLP_0068
do
    echo "===== $sample ====="
    samtools flagstat -@ "$SLURM_CPUS_PER_TASK" \
        "${sample}.sorted.bam" \
        > "${sample}.flagstat.txt"
done

echo "QC complete."
