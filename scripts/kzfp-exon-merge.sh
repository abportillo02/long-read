#!/bin/bash

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

bedtools="/home/abportillo/.conda/envs/mamba_abner_BC/bin/bedtools"

bedtools sort -i /home/abportillo/github_repo/long-read/scripts/kzfp_exons.bed | \
bedtools merge -i - -c 4 -o distinct > /home/abportillo/github_repo/long-read/docs/kzfp_exons_merged.bed
