#!/bin/bash
#SBATCH --job-name=plot_methylation
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=plot_methylation_%A_%a.out
#SBATCH --error=plot_methylation_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org


set -Eeuo pipefail

# Project locations
PROJECT="/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test"
SCRIPT="${PROJECT}/scripts/plot_methylation.py"
OUTPUT_DIR="${PROJECT}/mapping/DMR_analysis/figures"
LOG_DIR="${PROJECT}/logs"

# Use the Python environment that is active when this script is launched.
PYTHON="${PYTHON:-python}"

# Create output and log directories.
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# Check that the plotting script exists.
if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: Python script not found:"
    echo "  ${SCRIPT}"
    echo "Place plot_methylation.py in ${PROJECT}/scripts/"
    exit 1
fi

# Check required Python packages before starting the full analysis.
"${PYTHON}" -c "import numpy, pandas, matplotlib" || {
    echo "ERROR: numpy, pandas, and matplotlib must be installed in the active Python environment."
    echo "Activate the appropriate environment, or set PYTHON=/path/to/python."
    exit 1
}

# Timestamped log file.
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="${LOG_DIR}/plot_methylation_${TIMESTAMP}.log"

echo "Starting methylation plotting pipeline"
echo "Project: ${PROJECT}"
echo "Python:  $("${PYTHON}" --version 2>&1)"
echo "Script:  ${SCRIPT}"
echo "Log:     ${LOG_FILE}"
echo

# Run unbuffered so progress messages are written to the log immediately.
# Output is shown in the terminal and saved to the log simultaneously.
if "${PYTHON}" -u "${SCRIPT}" 2>&1 | tee "${LOG_FILE}"; then
    echo
    echo "SUCCESS: plotting pipeline finished."
    echo "Figures: ${OUTPUT_DIR}"
    echo "Log:     ${LOG_FILE}"
else
    status=$?
    echo
    echo "ERROR: plotting pipeline failed (exit status ${status})."
    echo "Check the log: ${LOG_FILE}"
    exit "${status}"
fi
