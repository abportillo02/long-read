#!/bin/bash
#SBATCH --job-name=FINALIZE_REFS
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=FINALIZE_%j.out
#SBATCH --error=FINALIZE_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -eo pipefail
source /home/abportillo/miniconda3/etc/profile.d/conda.sh
conda activate mamba_abner_BC
set -u

BASE=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference
GRCM39=$BASE/GRCm39.fa
cd $BASE/sample-refs
mkdir -p intermediate final scripts

# Junction strings (loxP is in reverse-complement orientation in the reads)
L5=TTCAATAACCCTTAAT
L3=CAGTCAGGTACATAAT
LOXRC=ATAACTTCGTATAATGTATGCTATACGAAGTTAT
R5=TAGGTCCCTCGAGG
R3=TAGGTGGATCCAGA

# 1. Move consensus outputs into intermediate/ (nothing is deleted)
for f in YRLP_00{66,67,68}.fa YRLP_00{66,67,68}.fa.fai YRLP_00{66,67,68}.chain YRLP_00{66,67,68}.consensus.log; do
  if [ -e "$f" ]; then mv "$f" intermediate/; fi
done

# 2. Save the patch script alongside the references
cat > scripts/patch_recomb.py << 'PY'
#!/usr/bin/env python3
import sys
L5="TTCAATAACCCTTAAT"; L3="CAGTCAGGTACATAAT"
LOXRC="ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
R5="TAGGTCCCTCGAGG"; R3="TAGGTGGATCCAGA"
fa_in, fa_out, chain_out = sys.argv[1:4]

def emit(header, seq, width, out):
    name = header[1:].split()[0]
    if name == "chr15":
        up = seq.upper()
        k5, k3 = L5+LOXRC+R5, L3+LOXRC+R3
        assert up.count(k5) == 1 and up.count(k3) == 1, "junctions not unique"
        p5 = up.index(k5) + len(L5) + len(LOXRC)   # cut starts at R5
        p3 = up.index(k3) + len(L3) + len(LOXRC)   # cut ends after 2nd loxP
        assert p5 < p3
        new = seq[:p5] + seq[p3:]
        with open(chain_out, "w") as c:
            c.write(f"chain 1000000 chr15 {len(seq)} + 0 {len(seq)} chr15 {len(new)} + 0 {len(new)} 1\n")
            c.write(f"{p5} {p3-p5} 0\n{len(seq)-p3}\n")
        print(f"cut chr15 {p5}-{p3} (0-based), deleted {p3-p5} bp, chr15 {len(seq)} -> {len(new)}")
        seq = new
    out.write(header + "\n")
    for i in range(0, len(seq), width):
        out.write(seq[i:i+width] + "\n")

with open(fa_in) as f, open(fa_out, "w") as out:
    header, buf, width = None, [], None
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                emit(header, "".join(buf), width, out)
            header, buf, width = line, [], None
        else:
            if width is None:
                width = len(line)
            buf.append(line)
    emit(header, "".join(buf), width, out)
PY

# 3. Patch 0066 chr15: floxed -> Cre-recombined (reads intermediate/, never edits it)
python3 scripts/patch_recomb.py \
  intermediate/YRLP_0066.fa \
  intermediate/YRLP_0066.recomb.fa \
  intermediate/YRLP_0066.floxed_to_recomb.chain

# 4. Final references, one clearly named file per sample
cp intermediate/YRLP_0067.fa        final/YRLP_0067_ref.fa
cp intermediate/YRLP_0068.fa        final/YRLP_0068_ref.fa
cp intermediate/YRLP_0066.recomb.fa final/YRLP_0066_ref.fa
for F in final/*_ref.fa; do samtools faidx $F; done

# 5. Checks (the job fails loudly if anything is off)
count() { samtools faidx "$1" chr15 | tail -n +2 | tr -d '\n' | tr a-z A-Z | { grep -F -o "$2" || true; } | wc -l; }

check() {  # file sample expected_floxed5 expected_floxed3 expected_recomb
  local a b c
  a=$(count "$1" "${L5}${LOXRC}${R5}")
  b=$(count "$1" "${L3}${LOXRC}${R3}")
  c=$(count "$1" "${L5}${LOXRC}${R3}")
  echo "$2 chr15 junctions: floxed5=$a floxed3=$b recombined=$c (expected $3 $4 $5)"
  if [ "$a" -ne "$3" ] || [ "$b" -ne "$4" ] || [ "$c" -ne "$5" ]; then echo "FAILED: $2"; exit 1; fi
}

check final/YRLP_0066_ref.fa YRLP_0066 0 0 1   # recombined
check final/YRLP_0067_ref.fa YRLP_0067 1 1 0   # floxed
check final/YRLP_0068_ref.fa YRLP_0068 0 0 0   # plain GRCm39

NREF=$(wc -l < ${GRCM39}.fai)
for S in YRLP_0066 YRLP_0067 YRLP_0068; do
  N=$(wc -l < final/${S}_ref.fa.fai)
  echo "$S contigs: $N (GRCm39: $NREF)"
  if [ "$N" -ne "$NREF" ]; then echo "FAILED: contig count $S"; exit 1; fi
done

echo "chr15 length, floxed -> recombined (0066):"
awk '$1=="chr15"{print "  before:", $2}' intermediate/YRLP_0066.fa.fai
awk '$1=="chr15"{print "  after: ", $2}' final/YRLP_0066_ref.fa.fai

# 6. Checksums and provenance
md5sum final/*_ref.fa intermediate/*.fa intermediate/*.chain > checksums.md5
cat > README.txt << EOF
Sample-specific references (generated $(date))
final/*_ref.fa   : references to use for mapping
  YRLP_0066 = GRCm39 + homozygous-ALT variants, then chr15 patched floxed -> Cre-recombined (one loxP, exon 3 interval removed)
  YRLP_0067 = GRCm39 + homozygous-ALT variants (floxed allele present as called)
  YRLP_0068 = GRCm39 + homozygous-ALT variants (CT26)
intermediate/    : bcftools consensus outputs, chain files (GRCm39 -> sample), floxed_to_recomb chain for 0066
scripts/         : patch_recomb.py
Variant filtering: GATK 4.1.8.0 hard filters, GQ>=20, DP floor/ceiling per sample (10/55, 12/70, 14/80), bcftools 1.24 norm, GT=1/1 only
Tools: $(samtools --version | head -n 1)
EOF

echo "DONE"