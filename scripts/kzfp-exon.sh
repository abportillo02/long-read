#!/bin/bash

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

gtf=$1 # /home/abportillo/github_repo/long-read/docs/gencode.v43.annotation.gtf
kzfp_genes=$2 # /home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv

grep -w "exon" "$gtf" | grep -F -f "$kzfp_genes" \
| awk 'BEGIN{OFS="\t"} {
    match($0, /gene_name "([^"]+)"/, a);
    print $1, $4-1, $5, a[1], ".", $7
}' > kzfp_exons.bed

