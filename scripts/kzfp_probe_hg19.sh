
#!/bin/bash
#SBATCH --job-name=probe_pipeline_hg19
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline_hg19.out
#SBATCH --error=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline_hg19.err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

# Activate conda environment
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

### INPUTS
GTF="/home/abportillo/github_repo/long-read/docs/gencode.v19.annotation.gtf"
FASTA="/home/abportillo/github_repo/long-read/docs/hg19.fa"
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"

### OUTPUT DIR
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila/hg19"
mkdir -p "$OUTDIR"

echo "STEP 1: Extract all exons for KZFP genes..."
awk '
BEGIN {
    while ((getline < "'"$KZFP_LIST"'") > 0) genes[$1]=1
}
$3=="exon" {
    match($0, /gene_name "([^"]+)"/, g)
    if (genes[g[1]]) print $0
}' "$GTF" > "$OUTDIR/kzfp_all_exons_hg19.gtf"
echo "STEP 1 complete"

echo "STEP 2: Convert exon GTF to BED..."
awk 'BEGIN {OFS="\t"} {
    if (match($0, /gene_name "([^"]+)"/, g)) gene=g[1]; else gene="NA"
    print $1, $4-1, $5, gene, ".", $7
}' "$OUTDIR/kzfp_all_exons_hg19.gtf" | sort -k1,1 -k2,2n > "$OUTDIR/kzfp_exons_raw_hg19.bed"
echo "STEP 2 complete"

echo "STEP 3: Assign exon numbers per gene..."
awk 'BEGIN {OFS="\t"} {
    if ($4 != prev_gene) { exon=1 } else { exon++ }
    print $1, $2, $3, $4"_exon"exon, ".", $6
    prev_gene=$4
}' "$OUTDIR/kzfp_exons_raw_hg19.bed" > "$OUTDIR/kzfp_exons_numbered_hg19.bed"
echo "STEP 3 complete"

echo "STEP 4: Extract exon sequences..."
bedtools getfasta -fi "$FASTA" -bed "$OUTDIR/kzfp_exons_numbered_hg19.bed" -name -s > "$OUTDIR/kzfp_exons_hg19.fa"
echo "STEP 4 complete"

echo "ALL DONE! Files in $OUTDIR"
