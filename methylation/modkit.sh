#!/bin/bash
#SBATCH --job-name=modkit-basic
#SBATCH --mail-type=ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=124G
#SBATCH --time=200:00:00
#SBATCH --mail-user=abportillo@coh.org
module load modkit

modkit pileup /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/NHA_Primary_p3_BR_Whole_C1_OUL14_L82343_PAQ01477_NNNNNNNN_L001.simplex.u.aligned.bam \ basic_mod_sample_PAQ01477.bedMethyl \
--log-filepath basic_pileup.log

module purge
