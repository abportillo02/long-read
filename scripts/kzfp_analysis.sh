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
REFERENCE="/home/abportillo/genomes/hg38/hg38_p14.fa"
ANNOTATION="/home/abportillo/genomes/hg38/gencode.v43.annotation.gtf"
BASE_OUTPUT_DIR="/home/abportillo/github_repo/long-read/results"

# Pool ALL fastq files under data/ into one file (outside FASTQ_DIR so it
# isn't picked up on re-runs). Handles both .fastq and .fastq.gz.
POOLED_DIR="/home/abportillo/github_repo/long-read/pooled_fastq"
POOLED="${POOLED_DIR}/MCF10A_pooled_pass.fastq"
mkdir -p "$POOLED_DIR"
> "$POOLED"   # start empty

find "$FASTQ_DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read -r f; do
    case "$f" in
        *.gz) zcat "$f" >> "$POOLED" ;;
        *)    cat  "$f" >> "$POOLED" ;;
    esac
done

# Sanity check: how many reads did we pool?
n_reads=$(( $(wc -l < "$POOLED") / 4 ))
echo "Pooled $n_reads reads into $POOLED"

# Run the pipeline once on the pooled sample
output_dir="${BASE_OUTPUT_DIR}/MCF10A_pooled"
echo "Processing pooled sample -> $output_dir"

python /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/dschones/KZFP_ISOFORM_ANALYSIS/KZFP_ISOFORM_PIPELINE/main.py \
    --fastq "$POOLED" \
    --reference "$REFERENCE" \
    --annotation "$ANNOTATION" \
    --output-dir "$output_dir"

echo "Done."