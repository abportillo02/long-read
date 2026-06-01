#!/bin/sh
#SBATCH --job-name=dorado_basecalling
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 1 
#SBATCH -N 1-4
#SBATCH -p gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=dorado_basecalling_%j.out
#SBATCH --error=dorado_basecalling_%j.err

# Load the necessary modules
module load Dorado/0.7.1

dorado basecaller \ --recursive \ --emit-moves\ --kit-name SQK-PCB114-24 \ --reference /scratch/abportillo/long-read/data/docs/hg38_p14.fa \ dna_r10,4.1_e8.2_400bps_sup@v5.0.0 \
/scratch/abportillo/long-read/data/data \ > /scratch/abportillo/long-read/data/aligned.bam 