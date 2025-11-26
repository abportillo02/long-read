#!/bin/bash
#SBATCH --job-name=probe_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline.out
#SBATCH --error=/home/abportillo/github_repo/long-read/kzfp_tequila/pipeline.err
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC
set -euo pipefail

### INPUTS 
GTF="/home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf"
FASTA="/home/abportillo/github_repo/long-read/docs/hg38_p14.fa"
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"   # one gene symbol per line

### OUTPUT DIR 
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila"
mkdir -p $OUTDIR

echo "STEP 1: Extract canonical transcript IDs..."

# Canonical transcript prioritization:
# MANE_Select > appris_principal > Ensembl_canonical
grep -F -f $KZFP_LIST $GTF | \
awk '
  $3=="transcript" && (
      $0 ~ /MANE_Select/ || 
      $0 ~ /appris_principal/ || 
      $0 ~ /Ensembl_canonical/
  ) {
      match($0, /transcript_id "([^"]+)"/, a)
      print a[1]
}' \
| sort -u \
> $OUTDIR/kzfp_canonical_transcripts.txt

echo "STEP 1 ✓  Canonical transcript list written."

echo "STEP 2: Extract canonical exons only..."

grep -F -f $OUTDIR/kzfp_canonical_transcripts.txt $GTF | \
awk '$3=="exon"' \
> $OUTDIR/kzfp_canonical_exons.gtf

echo "STEP 2 ✓  Extracted canonical exons."

echo "STEP 3: Convert exon GTF → BED..."

awk '
BEGIN {OFS="\t"}
{
    match($0, /gene_name "([^"]+)"/, g)
    gene=g[1]
    print $1, $4-1, $5, gene, ".", $7
}' \
$OUTDIR/kzfp_canonical_exons.gtf \
| sort -k1,1 -k2,2n -k3,3n \
> $OUTDIR/kzfp_exons_raw.bed

echo "STEP 3 ✓  BED file created."

echo "STEP 4: Assign exon numbers per gene..."

awk '
BEGIN {OFS="\t"}
{
  if ($4 != prev_gene) { exon = 1 }
  else { exon++ }

  print $1, $2, $3, $4"_exon"exon, ".", $6
  prev_gene = $4
}' \
$OUTDIR/kzfp_exons_raw.bed \
> $OUTDIR/kzfp_exons_numbered.bed

echo "STEP 4 ✓  Exon numbering complete."

echo "STEP 5: Extract exon sequences to FASTA..."

bedtools getfasta \
  -fi $FASTA \
  -bed $OUTDIR/kzfp_exons_numbered.bed \
  -name \
> $OUTDIR/kzfp_exons.fa

echo "STEP 5 ✓  FASTA extraction complete."

echo "ALL DONE!"
echo "Output directory: $OUTDIR"
echo "Files generated:"
echo " - kzfp_canonical_transcripts.txt"
echo " - kzfp_canonical_exons.gtf"
echo " - kzfp_exons_raw.bed"
echo " - kzfp_exons_numbered.bed"
echo " - kzfp_exons.fa"
