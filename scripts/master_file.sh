#!/bin/bash

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.12.01

# This script generates gene-based, master-transcript-based, transcript-based, or exon-based FASTA files from a GTF annotation file and a genome FASTA file.

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC


# Usage check
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <annotation_file> <genome_file> <gene_ID_list> <output_FASTA> <feature_type>"
    echo "Feature type must be one of: gene, transcript, exon, master"
    exit 1
fi

# Assign input arguments
annotation_file="$1"
genome_file="$2"
gene_id_list="$3"
output_fasta="$4"
feature_type="$5"

# Validate feature type
if [[ "$feature_type" != "gene" && "$feature_type" != "transcript" && "$feature_type" != "exon" && "$feature_type" != "master" ]]; then
    echo "Error: Invalid feature type '$feature_type'. Choose from: gene, transcript, exon, master"
    exit 1
fi

# Create a temporary directory
temp_dir=$(mktemp -d)
echo "Using temporary directory: $temp_dir"


# If feature type is "gene"
if [[ "$feature_type" == "gene" ]]; then
    echo "Filtering for gene features."

    # Build a temporary pattern file for gene_name attributes
    awk '{print "gene_name \""$0"\""}' "$gene_id_list" > "$temp_dir/gene_name_patterns.txt"

    # Filter the GTF file for lines with specified gene names and feature type "gene"
    grep -Ff "$temp_dir/gene_name_patterns.txt" "$annotation_file" | awk '$3 == "gene"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"
    
    # Convert filtered GTF data to BED
    gtf2bed < "$temp_dir/filtered.gtf" > "$temp_dir/temp.bed"
    bedtools sort -i "$temp_dir/temp.bed" > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    # Extract sequences
    bedtools getfasta -fi "$genome_file" -bed "$temp_dir/temp.sorted.bed" -name -s -fo "$output_fasta"
    echo "Sequences have been extracted and saved to $output_fasta."

# If feature type is "master"
elif [[ "$feature_type" == "master" ]]; then
    echo "Filtering for master transcript features."

    # Filter the GTF file for lines with specified gene IDs and feature type "exon"
    grep -Ff "$gene_id_list" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"

    # Convert filtered GTF data to BED
    gtf2bed < "$temp_dir/filtered.gtf" > "$temp_dir/temp.bed"
    bedtools sort -i "$temp_dir/temp.bed" > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    # Merge overlapping exons with gene ID information
    bedtools merge -i "$temp_dir/temp.sorted.bed" -c 4 -o distinct | sort -k1,1 -k2,2n > "$temp_dir/temp.merged.bed"
    echo "Overlapping exons merged."

    # Iterate through unique geneIDs in the merged bed file and extract sequences
    temp_gene_bed="$temp_dir/temp.merged.singlegene.bed"
            # for each unique gene ID in the merged bed file
            awk '{print $4}' "$temp_dir/temp.merged.bed" | sort -u | while read -r geneID; do
            # extract all regions for the current gene ID into the temporary bed file
            awk -v gene="$geneID" '$4 == gene' "$temp_dir/temp.merged.bed" > "$temp_gene_bed"
            # echo a single header line containing the gene ID to the output fasta file
            echo ">$geneID" >> "$output_fasta"
            # extract sequences for each interval without header lines or new line characters
            bedtools getfasta -fi "$genome_file" -bed "$temp_gene_bed" -fo | grep -v "^>" | tr -d '\n ' >> "$output_fasta"
            # move to a new line
            echo >> "$output_fasta"
            echo "Sequences obtained for $geneID"
    done
    echo "Sequences have been extracted and saved to $output_fasta."

# If feature type is "transcript"
elif [[ "$feature_type" == "transcript" ]]; then
    echo "Filtering for transcript features."

    # Filter the GTF file for lines with specified gene IDs and feature type "exon"
    grep -Ff "$gene_id_list" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"

    # Convert filtered GTF data to BED with transcript IDs
    awk 'BEGIN {OFS="\t"}
        $3 == "exon" {
            match($0, /transcript_id "([^"]+)"/, t)
            if (t[1] != "") {
                print $1, $4 - 1, $5, t[1], ".", $7
            }
        }' "$temp_dir/filtered.gtf" | sort -k1,1 -k2,2n > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    # Iterate through unique transcript ID in the bed file and extract sequences
    temp_transcript_bed="$temp_dir/temp.transcript.bed"
            # for each unique transcript ID in the bed file
            awk '{print $4}' "$temp_dir/temp.sorted.bed" | sort -u | while read -r transcriptID; do
            # extract all regions for the current transcript ID into the temporary bed file
            awk -v transcript="$transcriptID" '$4 == transcript' "$temp_dir/temp.sorted.bed" > "$temp_transcript_bed"
            # echo a single header line containing the transcript ID to the output fasta file
            echo ">$transcriptID" >> "$output_fasta"
            # extract sequences for each interval without header lines or new line characters
            bedtools getfasta -fi "$genome_file" -bed "$temp_transcript_bed" -fo | grep -v "^>" | tr -d '\n ' >> "$output_fasta"
            # move to a new line
            echo >> "$output_fasta"
            echo "Sequences obtained for $transcriptID"
    done
    echo "Sequences have been extracted and saved to $output_fasta."

# If feature type is "exon"
elif [[ "$feature_type" == "exon" ]]; then
    echo "Filtering for exon features."

    # Filter the GTF file for lines with specified gene IDs and feature type "exon"
    grep -Ff "$gene_id_list" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"

    # Convert filtered GTF data to BED
    gtf2bed < "$temp_dir/filtered.gtf" > "$temp_dir/temp.bed"
    bedtools sort -i "$temp_dir/temp.bed" > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    # Extract sequences and remove duplicates
    bedtools getfasta -fi "$genome_file" -bed "$temp_dir/temp.sorted.bed" -name -s -fo "$temp_dir/temp.sorted.fasta"
    seqkit rmdup -s "$temp_dir/temp.sorted.fasta" -o "$output_fasta"
    echo "Sequences have been extracted and saved to $output_fasta."

fi

# Remove the temporary files
rm -r "$temp_dir"

### Format to run Inputs (5 arguments):

# annotation_file — A GTF file (e.g., gene/transcript/exon annotations). Gnecode annotation file 
# genome_file — A genome FASTA (e.g., hg38.fa). Fasta file download 
# gene_ID_list — A plain text file, one gene ID per line (IDs must appear in the GTF lines). KZFP_list 
# output_fasta — Path to write extracted sequences. 
# feature_type — One of: gene, transcript, exon, master. either gene or master 


# bash your_script.sh annotation.gtf hg38.fa gene_ids.txt output.fasta gene

#bash master_file.sh /home/abportillo/genomes/hg38/gencode.v43.annotation.gtf /home/abportillo/genomes/hg38/hg38_p14.fa /home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv /home/abportillo/github_repo/long-read/kzfp_tequila/probes_hg38.fasta gene