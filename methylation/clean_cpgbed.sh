#!/bin/sh
module load htslib/1.17

#sort
sort -k1,1 -k2,2n /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.bed > /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.sorted.bed

#remove any old compressed versions
# rm -f cpgIslands-GRCh38.sorted.bed.gz*

#compress with bgzip (force overwrite)
bgzip -c /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.sorted.bed > /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.sorted.bed.gz

#index with tabix
tabix -p bed /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.sorted.bed.gz

#confirm it's BGZF
file /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/methylation_test/output/cpgIslands-GRCh38.sorted.bed.gz 


