#!/bin/bash

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

bedtools="/home/abportillo/.conda/envs/mamba_abner_BC/bin/bedtools"

bedtools getfasta \
  -fi /home/abportillo/github_repo/long-read/docs/hg38_p14.fa \
  -bed /home/abportillo/github_repo/long-read/scripts/kzfp_exons.bed \
  -name \
  -s \
  > /home/abportillo/github_repo/long-read/docs/kzfp_exons.fa
