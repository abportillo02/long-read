#!/bin/bash
#SBATCH --job-name=offtarget_intersect
#SBATCH -p all
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=offtarget_intersect%j.out
#SBATCH --error=offtarget_intersect%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

cd /home/abportillo/github_repo/kzfp_long_read/results/kzfp_signal/

# annotated, ranked, with a header — the keeper
{ echo -e "region\tmax_depth\tgene"
  bedtools intersect -a /home/abportillo/github_repo/kzfp_long_read/results/kzfp_signal/top_offtarget.bed -b /home/abportillo/github_repo/kzfp_long_read/results/kzfp_signal/gencode.v43.genes.bed -wa -wb \
    | sort -k4,4nr \
    | awk '{print $1":"$2"-"$3"\t"$4"\t"$NF}'
} > offtarget_annotated.tsv

# peek
# column -t offtarget_annotated.tsv | head -30