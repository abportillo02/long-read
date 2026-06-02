#!/bin/bash
# =============================================================================
# KZFP pilot diagnostics -- distinguishes "capture did nothing" from a script
# artifact hiding real signal. Run on apollo after the alignment step.
#
# Usage:  bash kzfp_diagnostics.sh
#         (activate the env first if samtools isn't on PATH:
#          conda activate /home/abportillo/.conda/envs/mamba_abner_BC )
#
# Writes everything to DIAGNOSTICS.txt in the results dir AND to the screen.
# =============================================================================

set -uo pipefail

DIR=/home/abportillo/github_repo/long-read/results/pilot_kzfp
cd "$DIR" || { echo "ERROR: cannot cd to $DIR" >&2; exit 1; }

OUT="$DIR/DIAGNOSTICS.txt"

# sanity: tools + inputs
if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not on PATH. Activate mamba_abner_BC first." >&2
    exit 1
fi
for f in panel.sorted.bam genome.sorted.bam pooled.fastq; do
    [[ -s "$DIR/$f" ]] || { echo "ERROR: missing $DIR/$f" >&2; exit 1; }
done

{
echo "==================== KZFP PILOT DIAGNOSTICS ===================="
echo "Date: $(date)"
echo "Dir:  $DIR"
echo

echo "### 1. Real per-gene counts (PRIMARY alignments only, de-inflated) ###"
echo "    [count] gene  -- compare to the inflated idxstats numbers"
samtools view -F 0x904 panel.sorted.bam | cut -f3 | sort | uniq -c | sort -rn | head -20
echo

echo "### 2. Panel MAPQ distribution (0 = ambiguous/repeat-smear, 60 = unique) ###"
echo "    [count] MAPQ"
samtools view -F 0x904 panel.sorted.bam | awk '{print $5}' | sort -n | uniq -c
echo

echo "### 3. Genome BED sanity + chrom-naming match ###"
echo "panel.genes.bed line count: $(wc -l < panel.genes.bed)"
echo "first 3 BED lines:"
head -3 panel.genes.bed
echo "chrom names actually present in genome.sorted.bam (mapped):"
samtools idxstats genome.sorted.bam | awk '$3>0{print $1}' | head
echo

echo "### 4. Where do the on-panel reads land on the GENOME? ###"
echo "    [count] chrom  MAPQ  -- scattered/non-KZFP chroms => artifact"
samtools view -F 0x904 panel.sorted.bam | cut -f1 | sort -u > panel_ids.txt
echo "unique reads with a primary panel alignment: $(wc -l < panel_ids.txt)"
samtools view -F 0x900 genome.sorted.bam \
  | awk 'NR==FNR{ids[$1]=1;next} ($1 in ids){print $3"\t"$5}' panel_ids.txt - \
  | sort | uniq -c | sort -rn | head -20
echo

echo "### 4b. How many on-panel reads map to the genome at ALL? ###"
mapped_on_genome=$(samtools view -F 0x904 genome.sorted.bam \
  | awk 'NR==FNR{ids[$1]=1;next} ($1 in ids)' panel_ids.txt - | cut -f1 | sort -u | wc -l)
echo "on-panel reads with a primary genome alignment: $mapped_on_genome / $(wc -l < panel_ids.txt)"
echo "(big gap => many 'panel' reads don't match the human genome = artifact)"
echo

echo "### 5. Read-length distribution (primer-dimer check) ###"
echo "    [count] length(bp)  -- spike at tens of bp = dimers"
awk 'NR%4==2{print length($0)}' pooled.fastq | sort -n | uniq -c
echo
echo "    length summary:"
awk 'NR%4==2{L=length($0); n++; s+=L; if(L<min||!min)min=L; if(L>max)max=L}
     END{printf "    reads=%d  min=%d  max=%d  mean=%.0f\n", n, min, max, s/n}' pooled.fastq
echo
echo "    binned:"
awk 'NR%4==2{L=length($0);
       if(L<100)b="<100";
       else if(L<300)b="100-299";
       else if(L<500)b="300-499";
       else if(L<1000)b="500-999";
       else if(L<2000)b="1000-1999";
       else b=">=2000";
       c[b]++}
     END{for(k in c) printf "    %-10s %d\n", k, c[k]}' pooled.fastq \
  | sort -t' ' -k1
echo
echo "==============================================================="
} | tee "$OUT"

echo
echo "Wrote: $OUT"