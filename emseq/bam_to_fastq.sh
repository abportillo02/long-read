#!/bin/bash
#SBATCH --job-name=custom_ref_methyl
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail

source activate mamba_abner_BC

# ---- Paths ----
BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
UBAM_DIR=$BASE/bam
OUTPUT_DIR=$BASE/fastq
mkdir -p "$OUTPUT_DIR"

# ---- Sample list: array index maps to BAM ----
BAMS=(
#   YRLP_0066_1_CO_Whole_C1_NFPFG_A48660_237232LT4_ACAGTGAC_L007_R.bam
#   YRLP_0067_1_CO_Whole_C1_NFPFG_A48661_237232LT4_TGACCAAC_L007_R.bam
  YRLP_0068_1_CL_Whole_C1_NFPFG_A48662_237232LT4_CGATGTAC_L007_R.bam
)

BAM=$UBAM_DIR/${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" .bam)
OUTDIR=$OUTPUT_DIR/${SAMPLE}
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "$OUTDIR/chrom_fa" logs

echo "[$SAMPLE] Starting pipeline on $(date)"

# ---- 1. Build per-chromosome consensus, then concatenate ----
# Splitting by chromosome makes this tractable at WGS scale since
# samtools consensus is single-threaded per invocation.
CHROMS=$(samtools view -H "$BAM" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}')

for CHR in $CHROMS; do
  samtools consensus -f fasta --min-depth 3 -r "$CHR" "$BAM" > "$OUTDIR/chrom_fa/${CHR}.fa" &
  # throttle to avoid oversubscribing the node
  while [ "$(jobs -r | wc -l)" -ge "$THREADS" ]; do wait -n; done
done
wait

cat "$OUTDIR"/chrom_fa/*.fa > "$OUTDIR/${SAMPLE}_consensus.fa"
samtools faidx "$OUTDIR/${SAMPLE}_consensus.fa"
echo "[$SAMPLE] Consensus reference built: $OUTDIR/${SAMPLE}_consensus.fa"

# ---- 2. Build minimap2 index once ----
minimap2 -d "$OUTDIR/${SAMPLE}_consensus.mmi" -x map-ont "$OUTDIR/${SAMPLE}_consensus.fa"

# ---- 3. Extract reads, preserving modification tags ----
samtools fastq -T MM,ML "$BAM" > "$OUTDIR/${SAMPLE}_reads.fq"

# ---- 4. Realign to this sample's own reference, carrying tags through ----
minimap2 -ax map-ont -y -t "$THREADS" \
  "$OUTDIR/${SAMPLE}_consensus.mmi" "$OUTDIR/${SAMPLE}_reads.fq" \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE}_realigned.sorted.bam"
samtools index "$OUTDIR/${SAMPLE}_realigned.sorted.bam"

# ---- 5. Sanity check: did MM/ML survive the round trip? ----
if ! samtools view "$OUTDIR/${SAMPLE}_realigned.sorted.bam" | head -1 | grep -q "MM:Z:"; then
  echo "[$SAMPLE] ERROR: MM tag missing after realignment. Check steps 3-4 before proceeding." >&2
  exit 1
fi
echo "[$SAMPLE] Modification tags confirmed present after realignment."

# ---- 6. Methylation calling ----
modkit pileup "$OUTDIR/${SAMPLE}_realigned.sorted.bam" "$OUTDIR/${SAMPLE}_methylation.bed" \
  -r "$OUTDIR/${SAMPLE}_consensus.fa" --threads "$THREADS"

echo "[$SAMPLE] Done on $(date). Output: $OUTDIR/${SAMPLE}_methylation.bed"

