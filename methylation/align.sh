#!/bin/bash
#SBATCH --job-name=align.sh    # Job name
#SBATCH --mail-type=END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=abportillo@coh.org     # Where to send mail  
#SBATCH -n 8                           # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=64G                      # Amount of memory in GB
#SBATCH --time=250:00:00               # Time limit hrs:min:sec
#SBATCH --output=dorado_align%j.log   # Standard output and error log

module load Dorado/0.5.2
module load samtools

#define variables
index="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/NHA_ONT/ONTsimplex/GRCh38.mmi"  
input_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/data"
output_directory="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output"

mkdir -p "$output_directory"

bam=(
    "NHA_Primary_p1_BR_Whole_C1_OGL14_L66007_PAK78455_NoIndex_L001_fb2f0c71.simplex.5mCG_5hmCG.ubam"
    "NHA_Primary_p3_BR_Whole_C1_OUL14_L82343_PAQ01477_NNNNNNNN_L001.simplex.u.bam"
    "NHA_Primary_p3_BR_Whole_C1_OUL14_L82343_PAO99537_NNNNNNNN_L001.simplex.u.bam")

#loop through files and align each one
for bam in "${bam[@]}"; do
    input_bam="$input_directory/$bam"
    
    #extract base name and remove only the last file extension (.ubam or .u.bam) since they have different extensions
    base_name=$(basename "$bam" | sed -E 's/\.[^.]+$//')

    output_bam="$output_directory/${base_name}.aligned.bam"

    echo "Aligning $input_bam..."
    dorado aligner "$index" "$input_bam" | samtools sort -o "$output_bam"
done
