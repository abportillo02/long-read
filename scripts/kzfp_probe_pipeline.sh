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

### INPUTS ###
GTF="/home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf"
FASTA="/home/abportillo/github_repo/long-read/docs/hg38_p14.fa"
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"   # one gene symbol per line

### OUTPUT DIR ###
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila"
mkdir -p $OUTDIR

echo "Extracting exon entries from GTF..."


# STEP 1 Extract all exons for the KZFP gene list

awk '$3=="exon"' $GTF | \
grep -F -f $KZFP_LIST | \
awk 'BEGIN{OFS="\t"} {
    # Extract chromosome, start, end, strand
    # Convert GTF 1-based to BED 0-based
    match($0, /gene_name "([^"]+)"/, a)
    gene=a[1]
    print $1, $4-1, $5, gene, ".", $7
}' \
| sort -k1,1 -k2,2n \
> $OUTDIR/kzfp_exons_raw.bed


echo "Merging exons per gene..."

# STEP 2 Merge exons per gene (unique biological exons)

bedtools merge \
    -i $OUTDIR/kzfp_exons_raw.bed \
    -c 4 -o distinct \
> $OUTDIR/kzfp_exons_merged.bed

echo "Adding exon numbers..."

# STEP 3 Add exon numbers per gene

awk '
BEGIN {OFS="\t"}
{
  if ($4 != prev_gene) { exon = 1 }
  else { exon++ }
  print $1, $2, $3, $4"_exon"exon, ".", "+"
  prev_gene = $4
}' $OUTDIR/kzfp_exons_merged.bed \
> $OUTDIR/kzfp_exons_numbered.bed

echo "Extracting exon FASTA..."


# STEP 4 Extract exon sequences from genome

bedtools getfasta \
    -fi $FASTA \
    -bed $OUTDIR/kzfp_exons_numbered.bed \
    -name \
> $OUTDIR/kzfp_exons.fa
EOF

echo "DONE!"
echo "Output directory: $OUTDIR"
echo "Files produced:"
echo " - kzfp_exons_raw.bed"
echo " - kzfp_exons_merged.bed"
echo " - kzfp_exons_numbered.bed"
echo " - kzfp_exons.fa"
