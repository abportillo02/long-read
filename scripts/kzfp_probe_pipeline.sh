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
KZFP_LIST="/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv"

### OUTPUT DIR ###
OUTDIR="/home/abportillo/github_repo/long-read/kzfp_tequila"
mkdir -p $OUTDIR

echo "Extracting exon entries from GTF..."

# STEP 1: Extract exons for the KZFP genes
awk '$3=="exon"' $GTF | \
grep -F -f $KZFP_LIST | \
awk 'BEGIN{OFS="\t"} {
    match($0, /gene_name "([^"]+)"/, a)
    gene=a[1]
    print $1, $4-1, $5, gene, ".", $7
}' | sort -k4,4 -k1,1 -k2,2n \
> $OUTDIR/kzfp_exons_raw.bed

echo "Merging exons PER GENE..."

# STEP 2: Merge exons for each gene separately
# Prevents multiple KZFP names from appearing in one row
awk '{
  print > ("'"$OUTDIR"'/tmp_"$4".bed")
}' $OUTDIR/kzfp_exons_raw.bed

# Merge each gene independently
for f in $OUTDIR/tmp_*.bed; do
    gene=$(basename $f | sed 's/tmp_//' | sed 's/.bed//')

    bedtools merge \
        -i $f \
        | awk -v g=$gene 'BEGIN{OFS="\t"} {print $1, $2, $3, g}'
done | sort -k4,4 -k1,1 -k2,2n \
> $OUTDIR/kzfp_exons_merged.bed

rm $OUTDIR/tmp_*.bed

echo "Adding exon numbers..."

# STEP 3: Add exon numbers
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

# STEP 4: Extract sequences
bedtools getfasta \
    -fi $FASTA \
    -bed $OUTDIR/kzfp_exons_numbered.bed \
    -name \
> $OUTDIR/kzfp_exons.fa

echo "DONE!"
echo "Output directory: $OUTDIR"
