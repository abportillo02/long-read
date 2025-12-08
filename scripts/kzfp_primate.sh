#!/bin/bash 
#SBATCH --job-name=kzfp-primate
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/docs/primate_list.out
#SBATCH --error=/home/abportillo/github_repo/long-read/docs/primate_list .err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --mem=150G
#SBATCH --time=48:00:00


# Activate conda environment
source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

python kzfp_primate.py
