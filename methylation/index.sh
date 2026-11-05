#!/bin/bash
#SBATCH --job-name=index.sh
#SBATCH --mail-type=ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=36G
#SBATCH --time=5:00:00
#SBATCH --mail-user=abportillo@coh.org
module load samtools
samtools index /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/NHA_Primary_p3_BR_Whole_C1_OUL14_L82343_PAQ01477_NNNNNNNN_L001.simplex.u.aligned.bam /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/NHA_Primary_p3_BR_Whole_C1_OUL14_L82343_PAQ01477_NNNNNNNN_L001.simplex.u.aligned.bam.bai
samtools faidx /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/files/GRCh38.fa
