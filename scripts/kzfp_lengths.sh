#!/bin/sh
#SBATCH --job-name=kzfp-lengths
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -N 1-1
#SBATCH --ntasks=1
#SBATCH --mem=75GB
#SBATCH --time=48:00:00
#SBATCH --output=kzfp_lengths.log

cd /home/abportillo/github_repo/long-read/docs

# 1. exact gene_name patterns (closing quote prevents ZNF8 matching ZNF80/ZNF808)
tr -d '\r' < priority_kzfps_ordered.txt \
  | awk '{gsub(/^[ \t]+|[ \t]+$/,"",$0)} $0!=""{print "gene_name \"" $0 "\""}' \
  > gene_patterns.txt

# 2. subset the GTF to your panel, protein-coding transcripts only
grep -Ff gene_patterns.txt gencode.v43.annotation.gtf \
  | awk -F'\t' '/transcript_type "protein_coding"/' > kzfp_panel_pc.gtf

# 3. extract spliced cDNA sequences (this is your proper FASTA)
gffread -w kzfp_cdna.fa -g hg38_p14.fa kzfp_panel_pc.gtf

# 4. lengths, longest first
seqkit fx2tab -nl kzfp_cdna.fa | sort -k2,2nr > kzfp_cdna_lengths.tsv
head kzfp_cdna_lengths.tsv