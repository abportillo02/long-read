#!/bin/bash

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
annotation_file="$1" # /home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf
genome_file="$2" # /home/abportillo/github_repo/long-read/docs/hg38_p14.fa
gene_id_list="$3" # /home/abportillo/github_repo/long-read/docs/primate_kzfp_protein_coding_matches_exact.txt
output_fasta="$4" # /home/abportillo/github_repo/long-read/docs/primate_kzfp_probes.fasta
feature_type="$5" # gene

# Validate feature type
if [[ "$feature_type" != "gene" && "$feature_type" != "transcript" && "$feature_type" != "exon" && "$feature_type" != "master" ]]; then
    echo "Error: Invalid feature type '$feature_type'. Choose from: gene, transcript, exon, master"
    exit 1
fi

# Create a temporary directory
temp_dir=$(mktemp -d)
echo "Using temporary directory: $temp_dir"

# Build gene_name patterns from CSV (column 1 = symbol), skip header, strip spaces/CRs
# If your list is plain text with NO header (one symbol per line), replace NR>1 with NR>=1 and drop sed.
awk -F',' 'NR>1 {sym=$1; gsub(/^[ \t]+|[ \t]+$/, "", sym); print "gene_name \""sym"\""}' "$gene_id_list" \
  | sed 's/\r$//' > "$temp_dir/gene_name_patterns.txt"


# # For TXT input: one symbol per line, trim spaces and CRLF
# tr -d '\r' < "$gene_id_list" \
#   | awk '
#       { sym=$0; gsub(/^[ \t]+|[ \t]+$/, "", sym) }
#       sym != "" { print "gene_name \"" sym "\""}
#     ' > "$temp_dir/gene_name_patterns.txt"


# If feature type is "gene"
if [[ "$feature_type" == "gene" ]]; then
    echo "Filtering for gene features."

    # Filter the GTF file for lines with specified gene names and feature type "gene"
    grep -Ff "$temp_dir/gene_name_patterns.txt" "$annotation_file" | awk '$3 == "gene"' > "$temp_dir/filtered.gtf"
    echo "Filtered rows (gene): $(wc -l < "$temp_dir/filtered.gtf")"

    if [[ $(wc -l < "$temp_dir/filtered.gtf") -eq 0 ]]; then
        echo "No 'gene' features matched your gene symbols. Check CSV column/header or symbols."
        rm -r "$temp_dir"; exit 2
    fi

    # Convert filtered GTF data to BED via gtf2bed (your original approach)
    gtf2bed < "$temp_dir/filtered.gtf" > "$temp_dir/temp.bed"

    # Map gene_id -> gene_name from the filtered GTF
    awk -F'\t' '
    $3 == "gene" {
        if (match($9, /gene_id "([^"]+)"/, gi) && match($9, /gene_name "([^"]+)"/, gn)) {
        print gi[1] "\t" gn[1]
        }
    }
    ' "$temp_dir/filtered.gtf" | sort -u > "$temp_dir/geneid_to_genename.tsv"

    # Replace BED column 4 (name) from gene_id -> gene_name
    awk 'BEGIN{FS=OFS="\t"}
        FNR==NR { m[$1]=$2; next }
        { if ($4 in m) $4=m[$4]; print }
    ' "$temp_dir/geneid_to_genename.tsv" "$temp_dir/temp.bed" > "$temp_dir/temp.named.bed"
 
    bedtools sort -i "$temp_dir/temp.named.bed" > "$temp_dir/temp.sorted.bed"
    echo "BED rows (gene): $(wc -l < "$temp_dir/temp.sorted.bed")"

    # Extract sequences
    bedtools getfasta -fi "$genome_file" -bed "$temp_dir/temp.sorted.bed" -name -s -fo "$output_fasta"
    echo "Sequences have been extracted and saved to $output_fasta."


elif [[ "$feature_type" == "master" ]]; then
    echo "Filtering for master transcript features."

    # Use gene_name patterns for exon filtering
    grep -Ff "$temp_dir/gene_name_patterns.txt" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered rows (exon): $(wc -l < "$temp_dir/filtered.gtf")"

    if [[ $(wc -l < "$temp_dir/filtered.gtf") -eq 0 ]]; then
        echo "No 'exon' features matched your gene symbols. Check CSV column/header or symbols."
        rm -r "$temp_dir"; exit 2
    fi

    # Build BED6 with gene_name in col4 and strand in col6
    awk 'BEGIN{OFS="\t"}
        $3=="exon"{
            match($0, /gene_name "([^"]+)"/, gn)
            if (gn[1] != "") print $1, $4 - 1, $5, gn[1], ".", $7
        }
    ' "$temp_dir/filtered.gtf" | bedtools sort -i - > "$temp_dir/temp.sorted.bed"
    echo "BED rows (exon BED for master): $(wc -l < "$temp_dir/temp.sorted.bed")"

    # >>> MINIMAL CHANGE: derive gene list from temp.sorted.bed (avoids combo labels like ZFP30,ZNF540)
    awk '{print $4}' "$temp_dir/temp.sorted.bed" | sort -u > "$temp_dir/gene_list.unique.txt"

    temp_gene_bed="$temp_dir/temp.gene.bed"
    : > "$output_fasta"

    # Iterate per gene; subset, merge within gene, and extract in transcription order
    while read -r geneID; do
        [[ -z "$geneID" ]] && continue

        # Subset only this gene from the pre-sorted BED6
        awk -v gene="$geneID" 'BEGIN{OFS="\t"} $4==gene {print $0}' "$temp_dir/temp.sorted.bed" > "$temp_gene_bed"
        [[ ! -s "$temp_gene_bed" ]] && echo "No exon intervals for $geneID" && continue

        # Capture strand (assumes consistent strand per gene; typical for protein-coding)
        strand=$(awk 'NR==1{print $6}' "$temp_gene_bed")

        # Merge within gene (strand-aware)
        bedtools merge -s -i "$temp_gene_bed" > "$temp_dir/${geneID}.merged.bed"

        # Order merged intervals by transcription direction
        if [[ "$strand" == "-" ]]; then
            sort -k1,1 -k2,2nr "$temp_dir/${geneID}.merged.bed" > "$temp_dir/${geneID}.ordered.bed3"
        else
            sort -k1,1 -k2,2n  "$temp_dir/${geneID}.merged.bed" > "$temp_dir/${geneID}.ordered.bed3"
        fi

        # Convert BED3 back to BED6 for strand-aware getfasta
        awk -v g="$geneID" -v s="$strand" 'BEGIN{OFS="\t"} {print $1,$2,$3,g,".",s}' \
            "$temp_dir/${geneID}.ordered.bed3" > "$temp_dir/${geneID}.ordered.bed6"

        # Write FASTA header and concatenated sequence (reverse-complement applied if '-' strand)
        echo ">$geneID" >> "$output_fasta"
        bedtools getfasta -fi "$genome_file" -bed "$temp_dir/${geneID}.ordered.bed6" -s \
            | grep -v "^>" | tr -d '\n ' >> "$output_fasta"
        echo >> "$output_fasta"

        echo "Sequences obtained for $geneID"
    done < "$temp_dir/gene_list.unique.txt"

    echo "Master sequences have been extracted and saved to $output_fasta."


# If feature type is "transcript"
elif [[ "$feature_type" == "transcript" ]]; then
    echo "Filtering for transcript features."

    grep -Ff "$gene_id_list" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"

    awk 'BEGIN {OFS="\t"}
        $3 == "exon" {
            match($0, /transcript_id "([^"]+)"/, t)
            if (t[1] != "") {
                print $1, $4 - 1, $5, t[1], ".", $7
            }
        }' "$temp_dir/filtered.gtf" | sort -k1,1 -k2,2n > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    temp_transcript_bed="$temp_dir/temp.transcript.bed"
    awk '{print $4}' "$temp_dir/temp.sorted.bed" | sort -u | while read -r transcriptID; do
        awk -v transcript="$transcriptID" '$4 == transcript' "$temp_dir/temp.sorted.bed" > "$temp_transcript_bed"
        echo ">$transcriptID" >> "$output_fasta"
        bedtools getfasta -fi "$genome_file" -bed "$temp_transcript_bed" | grep -v "^>" | tr -d '\n ' >> "$output_fasta"
        echo >> "$output_fasta"
        echo "Sequences obtained for $transcriptID"
    done
    echo "Sequences have been extracted and saved to $output_fasta."

# If feature type is "exon"
elif [[ "$feature_type" == "exon" ]]; then
    echo "Filtering for exon features."

    grep -Ff "$gene_id_list" "$annotation_file" | awk '$3 == "exon"' > "$temp_dir/filtered.gtf"
    echo "Filtered GTF data saved to $temp_dir/filtered.gtf"

    gtf2bed < "$temp_dir/filtered.gtf" > "$temp_dir/temp.bed"
    bedtools sort -i "$temp_dir/temp.bed" > "$temp_dir/temp.sorted.bed"
    echo "Sorted BED saved to $temp_dir/temp.sorted.bed"

    bedtools getfasta -fi "$genome_file" -bed "$temp_dir/temp.sorted.bed" -name -s -fo "$temp_dir/temp.sorted.fasta"
    seqkit rmdup -s "$temp_dir/temp.sorted.fasta" -o "$output_fasta"
    echo "Sequences have been extracted and saved to $output_fasta."
fi

# Remove the temporary files
rm -r "$temp_dir"
