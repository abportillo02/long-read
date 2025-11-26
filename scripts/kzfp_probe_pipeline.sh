
#!/bin/bash
#SBATCH --job-name=probe_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline.out
#SBATCH --error=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline.err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

set -euo pipefail

# Activate environment
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

### INPUTS
GTF="/home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf"
FASTA="/home/abportillo/github_repo/long-read/docs/hg38_p14.fa"
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"   # one gene symbol per line

### OUTPUT DIR
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila"
mkdir -p "$OUTDIR"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"; }

# Check inputs
for f in "$GTF" "$FASTA" "$KZFP_LIST"; do
    [[ -f "$f" ]] || { echo "ERROR: Missing file $f"; exit 1; }
done

log "STEP 1: Extract canonical transcript IDs..."
awk '
BEGIN {
    while ((getline < "'$KZFP_LIST'") > 0) genes[$1]=1
}
$3=="transcript" && (
    $0 ~ /MANE_Select/ ||
    $0 ~ /appris_principal/ ||
    $0 ~ /Ensembl_canonical/
) {
    match($0, /gene_name "([^"]+)"/, g)
    match($0, /transcript_id "([^"]+)"/, t)
    if (genes[g[1]]) print t[1]
}' "$GTF" | sort -u > "$OUTDIR/kzfp_canonical_transcripts.txt"
log "STEP 1 ✓ Canonical transcript list written."

log "STEP 2: Extract canonical exons only..."
awk '
BEGIN {
    while ((getline < "'$OUTDIR/kzfp_canonical_transcripts.txt'") > 0) tx[$1]=1
}
$3=="exon" {
    match($0, /transcript_id "([^"]+)"/, t)
    if (tx[t[1]]) print $0
}' "$GTF" > "$OUTDIR/kzfp_canonical_exons.gtf"
log "STEP 2 ✓ Extracted canonical exons."

log "STEP 3: Convert exon GTF → BED..."
awk '
BEGIN {OFS="\t"}
{
    if (match($0, /gene_name "([^"]+)"/, g)) gene=g[1]; else gene="NA"
    print $1, $4-1, $5, gene, ".", $7
}' "$OUTDIR/kzfp_canonical_exons.gtf" | sort -k1,1 -k2,2n > "$OUTDIR/kzfp_exons_raw.bed"
log "STEP 3 ✓ BED file created."

log "STEP 4: Assign exon numbers per gene..."
awk '
BEGIN {OFS="\t"}
{
    if ($4 != prev_gene) { exon = 1 } else { exon++ }
    print $1, $2, $3, $4"_exon"exon, ".", $6
    prev_gene = $4
}' "$OUTDIR/kzfp_exons_raw.bed" > "$OUTDIR/kzfp_exons_numbered.bed"
log "STEP 4 ✓ Exon numbering complete."

log "STEP 5: Extract exon sequences to FASTA..."
bedtools getfasta -fi "$FASTA" -bed "$OUTDIR/kzfp_exons_numbered.bed" -name -s > "$OUTDIR/kzfp_exons.fa"
log "STEP 5 ✓ FASTA extraction complete."

log "ALL DONE!"
log "Output directory: $OUTDIR"
log "Files generated:"
log " - kzfp_canonical_transcripts.txt"
log " - kzfp_canonical_exons.gtf"
log " - kzfp_exons_raw.bed"
log " - kzfp_exons_numbered.bed"
log " - kzfp_exons.fa"
