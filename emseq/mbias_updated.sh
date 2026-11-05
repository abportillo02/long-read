#!/bin/bash
#SBATCH --job-name=mbias
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=mbias_%j.out
#SBATCH --error=mbias_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=abportillo@coh.org

set -euo pipefail

export PATH="/home/abportillo/.local/share/mamba/envs/methyldackel_env/bin:$PATH"

EM="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test"
REFDIR="$EM/genotype-reference/sample-refs/mapping-refs"
TABLE="$EM/mapping/libraries.tsv"

cd "$EM/mapping"

# Preserve existing results rather than overwrite or delete them.
if [[ -d mbias ]]; then
    BACKUP="mbias_previous_$(date +%Y%m%d_%H%M%S)"
    mv mbias "$BACKUP"
    echo "Previous results preserved in $BACKUP"
fi

mkdir -p mbias

date
echo "MethylDackel: $(MethylDackel --version 2>&1 | head -n 1 || true)"
echo "Libraries file: $TABLE"

[[ -s "$TABLE" ]] || {
    echo "ERROR: Missing or empty libraries.tsv: $TABLE"
    exit 1
}

# Preserve the original file layout:
# column 1 = library ID; column 4 = reference ID.
# Skip the header only if the first row is a header.
while IFS=$'\t' read -r L _ _ R _; do
    [[ -n "${L:-}" && -n "${R:-}" ]] || continue

    # Skip a header row if present.
    [[ "$L" == "library" || "$L" == "Library" ]] && continue

    BAM="${L}.markdup.bam"
    REF="$REFDIR/${R}_ref.withcontrols.fa"

    [[ -s "$BAM" ]] || {
        echo "ERROR: Missing BAM for $L: $BAM"
        exit 1
    }

    [[ -s "$REF" ]] || {
        echo "ERROR: Missing reference for $L: $REF"
        exit 1
    }

    echo "=========================================="
    echo "M-bias: $L (reference $R)"
    date

    # Genome-wide M-bias; no -r chr1 restriction.
    # --txt outputs tabular metrics to stdout.
    # SVG plots remain enabled.
    MethylDackel mbias \
        --txt \
        -@ "${SLURM_CPUS_PER_TASK:-2}" \
        "$REF" \
        "$BAM" \
        "mbias/${L}" \
        > "mbias/${L}.mbias.tsv" \
        2> "mbias/${L}.mbias.log"

    [[ -s "mbias/${L}.mbias.tsv" ]] || {
        echo "ERROR: Empty M-bias TSV for $L"
        exit 1
    }

    echo "Completed: $L"
done < "$TABLE"

date
echo "DONE"
