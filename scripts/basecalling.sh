#!/bin/bash
#SBATCH --job-name=kzfp_basecall
#SBATCH --partition=gpu-a100             
#SBATCH --gres=gpu:1                
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=kzfp_basecall_%j.out
#SBATCH --error=kzfp_basecall_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

# =============================================================================
# STEP 1 of 2 -- BASECALLING ONLY (GPU, run on gemini)
#
# All 75 pod5 (barcodes 02/03/unclassified, all MCF10A) are basecalled POOLED.
# Barcodes are trimmed but NOT demultiplexed. Output is one unaligned BAM.
# Dorado writes BAM natively, so this step needs only Dorado -- no samtools.
#
# Hand off calls.bam to the apollo alignment script (step 2).
# =============================================================================

set -uo pipefail

# --- Environment ---
module load Dorado/0.7.1    


# --- Paths (EDIT) ---
BASE=/scratch/abportillo/long-read
POD5_DIR=$BASE/data/pod5                           
MODEL=MODEL="/scratch/abportillo/long-read/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
OUTDIR=$BASE/results/pilot_kzfp
mkdir -p "$OUTDIR"

echo "[$(date)] Basecalling on $(hostname)..."
${DORADO:-dorado} basecaller \
    -x cuda:all \
    --recursive \
    --kit-name SQK-PCB114-24 \
    --min-qscore 3 \
    "$MODEL" \
    "$POD5_DIR" \
  > "$OUTDIR/calls.bam"

n=$(${DORADO:-dorado} summary "$OUTDIR/calls.bam" 2>/dev/null | tail -n +2 | wc -l || echo "?")
echo "[$(date)] Done. Wrote $OUTDIR/calls.bam (approx $n read records)"
echo
echo "NEXT: run 2_kzfp_align_apollo.slurm.sh on apollo against:"
echo "  $OUTDIR/calls.bam"
echo "If gemini and apollo do NOT share this filesystem, copy calls.bam to apollo first."