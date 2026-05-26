#!/bin/bash
#SBATCH --job-name=kzfp_analysis
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/docs/kzfp_analysis.out
#SBATCH --error=/home/abportillo/github_repo/long-read/docs/kzfp_analysis.err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

# Configuration
FASTQ_DIR="/home/abportillo/github_repo/long-read/data"
FASTQ_LIST="/home/abportillo/github_repo/long-read/scripts/fastq_list.txt"
REFERENCE="/home/abportillo/genomes/hg38/hg38_p14.fa"
ANNOTATION="/home/abportillo/genomes/hg38/gencode.v43.annotation.gtf"
BASE_OUTPUT_DIR="/home/abportillo/github_repo/long-read/results"

# Loop through each fastq file in the list
while IFS= read -r fastq_file; do
    # Skip empty lines
    [[ -z "$fastq_file" ]] && continue

    # Create a unique output directory per sample (strip .fastq extension)
    sample_name="${fastq_file%.fastq}"
    output_dir="${BASE_OUTPUT_DIR}/${sample_name}"

    echo "Processing: $fastq_file -> $output_dir"

    python /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/dschones/KZFP_ISOFORM_ANALYSIS/KZFP_ISOFORM_PIPELINE/main.py \
        --fastq "${FASTQ_DIR}/${fastq_file}" \
        --reference "$REFERENCE" \
        --annotation "$ANNOTATION" \
        --output-dir "$output_dir"

done < "$FASTQ_LIST"

echo "All samples processed."