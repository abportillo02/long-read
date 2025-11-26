#!/bin/bash

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

input="/home/abportillo/github_repo/long-read/scripts/kzfp_exons.bed"
output="/home/abportillo/github_repo/long-read/docs/kzfp_exons_numbered.bed"

awk '{
    key=$4;
    if(key!=prev){count=1} else {count++}
    prev=key;
    print $1, $2, $3, key"_exon"count, $5
}' OFS="\t" <(sort -k4,4 -k2,2n $input) > $output
