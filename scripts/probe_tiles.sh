#!/bin/bash 
#SBATCH --job-name=probe_tiling
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH --output=/home/abportillo/github_repo/long-read/kzfp_tequila/probe_tiling.out
#SBATCH --error=/home/abportillo/github_repo/long-read/kzfp_tequila/probe_tiling.err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --mem=150G
#SBATCH --time=48:00:00

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

## Probe tiling python script

python probe_tiles.py /home/abportillo/github_repo/long-read/kzfp_tequila/kzfp_exons.fa /home/abportillo/github_repo/long-read/kzfp_tequila/kzfp_tequila_probes.csv