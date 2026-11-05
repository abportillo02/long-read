for f in variant-calling/vcf/*.raw.vcf.gz
do
    echo "===== $f ====="
    gatk VariantsToTable \
        -V "$f" \
        -F CHROM \
        -F POS \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -GF GT \
        -GF DP \
        -GF GQ \
        -O "${f%.vcf.gz}.summary.tsv"
done