#!/bin/bash
#SBATCH --job-name=methyl_compare
#SBATCH --array=1-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=36G
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/methyl_compare_%A_%a.log
#SBATCH --error=logs/methyl_compare_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

