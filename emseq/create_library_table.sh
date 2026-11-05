EM=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test
mkdir -p $EM/mapping

printf '%s\t%s\t%s\t%s\n' \
 YRLP_0064 5hmC-Seq AMPK-KO-CD  YRLP_0066 \
 YRLP_0065 EM-Seq   AMPK-WT-HFD YRLP_0067 \
 YRLP_0069 EM-Seq   AMPK-KO-HFD YRLP_0066 \
 YRLP_0070 EM-Seq   CT26        YRLP_0068 \
 YRLP_0071 5hmC-Seq AMPK-WT-HFD YRLP_0067 \
 YRLP_0072 5hmC-Seq AMPK-KO-HFD YRLP_0066 \
 YRLP_0073 5hmC-Seq CT26        YRLP_0068 \
 YRLP_0074 5hmC-Seq AMPK-WT-CD  YRLP_0067 \
 YRLP_0075 EM-Seq   AMPK-WT-CD  YRLP_0067 \
 YRLP_0076 EM-Seq   AMPK-KO-CD  YRLP_0066 \
 > $EM/mapping/libraries.tsv

ls -lh $EM/mapping/libraries.tsv
awk -F'\t' '{print NF}' $EM/mapping/libraries.tsv | sort -u     # should print only 4
cat $EM/mapping/libraries.tsv

FQ=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/20251029_LH00295_0299_B237232LT4/raw-data/fastq
for L in $(cut -f1 $EM/mapping/libraries.tsv); do
  echo "$L R1: $(ls $FQ/${L}_*_R1_001.fastq.gz | wc -l)  R2: $(ls $FQ/${L}_*_R2_001.fastq.gz | wc -l)"
done