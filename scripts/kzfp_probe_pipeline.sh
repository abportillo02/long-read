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

# Activate conda environment (no bashrc sourcing)
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

### INPUTS
GTF="/home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf"
FASTA="/home/abportillo/github_repo/long-read/docs/hg38_p14.fa"
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"

### OUTPUT DIR
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila"
mkdir -p "$OUTDIR"

echo "STEP 1: Extract canonical transcript IDs..."
awk '
BEGIN {
    while ((getline < "'"$KZFP_LIST"'") > 0) genes[$1]=1
}
$3=="transcript" && ($0 ~ /MANE_Select/) {
    match($0, /gene_name "([^"]+)"/, g)
    match($0, /transcript_id "([^"]+)"/, t)
    if (genes[g[1]]) print t[1]
}' "$GTF" | sort -u > "$OUTDIR/kzfp_canonical_transcripts.txt"
echo "STEP 1 ✓"

echo "STEP 2: Extract canonical exons..."
awk '
BEGIN {
    while ((getline < "'"$OUTDIR/kzfp_canonical_transcripts.txt"'") > 0) tx[$1]=1
}
$3=="exon" {
    match($0, /transcript_id "([^"]+)"/, t)
    if (tx[t[1]]) print $0
}' "$GTF" > "$OUTDIR/kzfp_canonical_exons.gtf"
echo "STEP 2 ✓"

echo "STEP 3: Convert exon GTF → BED..."
awk 'BEGIN {OFS="\t"} {
    if (match($0, /gene_name "([^"]+)"/, g)) gene=g[1]; else gene="NA"
    print $1, $4-1, $5, gene, ".", $7
}' "$OUTDIR/kzfp_canonical_exons.gtf" | sort -k1,1 -k2,2n > "$OUTDIR/kzfp_exons_raw.bed"
echo "STEP 3 ✓"

echo "STEP 4: Assign exon numbers..."
awk 'BEGIN {OFS="\t"} {
    if ($4 != prev_gene) { exon=1 } else { exon++ }
    print $1, $2, $3, $4"_exon"exon, ".", $6
    prev_gene=$4
}' "$OUTDIR/kzfp_exons_raw.bed" > "$OUTDIR/kzfp_exons_numbered.bed"
echo "STEP 4 ✓"

echo "STEP 5: Extract exon sequences..."
bedtools getfasta -fi "$FASTA" -bed "$OUTDIR/kzfp_exons_numbered.bed" -name -s > "$OUTDIR/kzfp_exons.fa"
echo "STEP 5 ✓"

echo "ALL DONE! Files in $OUTDIR"
