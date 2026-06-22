#!/bin/sh
#SBATCH --job-name=transfer_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -N 1-1
#SBATCH --ntasks=1
#SBATCH --mem=75GB
#SBATCH --time=48:00:00
#SBATCH --output=transfer_%A_%a.log

scp -r /home/abportillo/github_repo/kzfp_long_read/data/pod5/ gemini-data1:/scratch/abportillo/long-read/data