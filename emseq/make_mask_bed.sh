#!/bin/bash
#SBATCH --job-name=maskbed
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/maskbed_%A_%a.out
#SBATCH --error=logs/maskbed_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

export PATH=/home/abportillo/.local/share/mamba/envs/bcftools_env/bin:$PATH
cd /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/variant-calling/filtered

for S in YRLP_0066 YRLP_0067 YRLP_0068; do
  bcftools view -H -v snps -i 'GT="alt" && ((REF="C" && ALT="T") || (REF="G" && ALT="A"))' \
    ${S}.pass.norm.vcf.gz \
    | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$4">"$5}' > ${S}.CtoT_GtoA.mask.bed
  wc -l ${S}.CtoT_GtoA.mask.bed
done
make_mask_bed.sh (END)

