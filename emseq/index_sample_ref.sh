## Ran on bash 
## After making mask bed files and sample references fasta files
## Sanity check 

cd /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/sample-refs

# index every FASTA (the reference may already have an index, -f is not needed for faidx)
for S in YRLP_0066 YRLP_0067 YRLP_0068; do
  samtools faidx ${S}.fa
done
samtools faidx ../GRCm39.fa

# contig counts and total length
for F in YRLP_0066.fa YRLP_0067.fa YRLP_0068.fa ../GRCm39.fa; do
  echo "== $F"
  awk '{s+=$2} END {print "contigs:", NR, " total bp:", s}' ${F}.fai
done

# contig names and lengths should match GRCm39 except for indel-driven length changes
for S in YRLP_0066 YRLP_0067 YRLP_0068; do
  echo "== $S: contig names identical to GRCm39?"
  diff <(cut -f1 ${S}.fa.fai) <(cut -f1 ../GRCm39.fa.fai) && echo "yes"
done

# chain files
ls -lh *.chain
head -n 3 YRLP_0066.chain