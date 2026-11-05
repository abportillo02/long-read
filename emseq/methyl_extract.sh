#!/bin/bash
#SBATCH --job-name=mdextract
#SBATCH --array=1-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=36G
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/mdextract_%A_%a.log
#SBATCH --error=logs/mdextract_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail

export PATH="/home/abportillo/.local/share/mamba/envs/methyldackel_env/bin:$PATH"

command -v MethylDackel || {
    echo "ERROR: MethylDackel not found in PATH" >&2
    exit 1
}

MethylDackel --version

MAP=/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/mapping
OUT=$MAP/methyl
LIBS=$MAP/libraries.tsv

REF_TEMPLATE="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/genotype-reference/sample-refs/mapping-refs/{REF}_ref.withcontrols.fa"

# Read the library assigned to this array task
IFS=$'\t' read -r LIB ASSAY GROUP REFID < <(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIBS"
)

[[ -n "$LIB" && -n "$ASSAY" && -n "$REFID" ]] || {
    echo "Missing library information on row $SLURM_ARRAY_TASK_ID" >&2
    exit 1
}

# Assay-specific M-bias trimming windows
case "$ASSAY" in
    EM-Seq)
        WIN=(--OT 3,150,0,129 --OB 2,149,23,0)
        ;;
    5hmC-Seq)
        WIN=(--OT 4,138,32,146 --OB 14,148,6,120)
        ;;
    *)
        echo "Unknown assay label: $ASSAY" >&2
        exit 1
        ;;
esac

REF=${REF_TEMPLATE/\{REF\}/$REFID}
BAM=$MAP/$LIB.markdup.bam

[[ -s "$REF" ]] || {
    echo "Reference not found or empty: $REF" >&2
    exit 1
}

[[ -s "$BAM" ]] || {
    echo "BAM not found or empty: $BAM" >&2
    exit 1
}

[[ -s "$BAM.bai" || -s "${BAM%.bam}.bai" ]] || {
    echo "BAM index not found for: $BAM" >&2
    exit 1
}

mkdir -p "$OUT"

{
    echo "date: $(date -Is)"
    echo "library: $LIB"
    echo "assay: $ASSAY"
    echo "group: $GROUP"
    echo "reference: $REF"
    echo "MethylDackel: $(MethylDackel --version 2>&1 | head -1)"
    printf 'M-bias options:'
    printf ' %s' "${WIN[@]}"
    echo
    echo "other options: --mergeContext --minDepth 1"
} > "$OUT/$LIB.extract.params.txt"

MethylDackel extract \
    "${WIN[@]}" \
    --mergeContext \
    --minDepth 1 \
    -@ "$SLURM_CPUS_PER_TASK" \
    -o "$OUT/$LIB" \
    "$REF" \
    "$BAM"

echo "Finished extraction for $LIB"
```
