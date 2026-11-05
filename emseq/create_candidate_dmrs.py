#!/usr/bin/env python3
"""
this is making exploratory candidate DMR tables for plotting from MethylDackel CpG bedGraphs.

Input files:
  /net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/mapping/methyl/

Each bedGraph contains:
  chromosome, start, end, methylation_percent, methylated_reads, unmethylated_reads

The script:
  1. Reads large bedGraph files in chunks (without assuming every row is valid).
  2. Removes malformed rows and CpGs below the coverage threshold.
  3. Summarizes CpGs into fixed genomic windows.
  4. Calculates the four planned contrasts and genotype-by-diet interaction.
  5. Writes all-window and candidate-window TSV files for EM-seq and 5hmC-seq.

These are exploratory candidate regions, not statistically validated DMRs.
There is one biological sample per condition, so biological variance and
inferential significance cannot be estimated.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

BASE = Path(
    "/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/mapping/methyl"
)
OUTDIR = BASE.parent / "DMR_analysis"

CHUNK_SIZE = 1_000_000
WINDOW_SIZE = 1_000          # 1 kb genomic windows
MIN_CPG_COVERAGE = 5         # minimum reads at an individual CpG
MIN_CPGS_PER_WINDOW = 5      # minimum retained CpGs in a window
MIN_ABS_DELTA = 10.0         # percentage-point difference to flag a candidate


ASSAYS = {
    "EM_seq": {
        "YRLP_0064": "KO_CD",
        "YRLP_0069": "KO_HFD",
        "YRLP_0075": "WT_CD",
        "YRLP_0065": "WT_HFD",
    },
    "5hmC_seq": {
        "YRLP_0076": "KO_CD",
        "YRLP_0072": "KO_HFD",
        "YRLP_0074": "WT_CD",
        "YRLP_0071": "WT_HFD",
    },
}

NUMERIC_COLUMNS = ["start", "end", "percent", "meth", "unmeth"]


# ---------------------------------------------------------------------
# Input processing
# ---------------------------------------------------------------------

def summarize_one_file(filepath: Path, sample_name: str) -> pd.DataFrame:
    """Read one large bedGraph and summarize retained CpGs by 1-kb window."""

    print(f"  Reading {filepath.name}", flush=True)

    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    # Read as strings first. This prevents pandas from crashing if a malformed
    # row contains a blank or non-numeric coordinate.
    reader = pd.read_csv(
        filepath,
        sep="\t",
        comment="t",       # skips MethylDackel track header
        header=None,
        names=["chrom", "start", "end", "percent", "meth", "unmeth"],
        chunksize=CHUNK_SIZE,
        dtype="string",
        on_bad_lines="skip",
        low_memory=False,
    )

    summaries = []
    rows_seen = 0
    rows_kept = 0

    for chunk in reader:
        rows_seen += len(chunk)

        # Convert numeric columns; malformed values become NaN.
        for col in NUMERIC_COLUMNS:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce")

        chunk["chrom"] = chunk["chrom"].astype("string").str.strip()

        # Remove blank/malformed records.
        chunk = chunk.dropna(
            subset=["chrom", "start", "end", "meth", "unmeth"]
        )
        chunk = chunk[chunk["chrom"] != ""]

        # Coordinates and read counts must be nonnegative integers.
        chunk = chunk[
            (chunk["start"] >= 0)
            & (chunk["end"] > chunk["start"])
            & (chunk["meth"] >= 0)
            & (chunk["unmeth"] >= 0)
        ].copy()

        chunk["start"] = chunk["start"].astype("int64")
        chunk["end"] = chunk["end"].astype("int64")
        chunk["meth"] = chunk["meth"].astype("int64")
        chunk["unmeth"] = chunk["unmeth"].astype("int64")

        chunk["coverage"] = chunk["meth"] + chunk["unmeth"]

        # Retain CpGs with sufficient read coverage.
        chunk = chunk[chunk["coverage"] >= MIN_CPG_COVERAGE].copy()
        if chunk.empty:
            continue

        rows_kept += len(chunk)

        # Use the CpG start coordinate to assign a fixed-width window.
        chunk["window_start"] = (
            chunk["start"] // WINDOW_SIZE
        ) * WINDOW_SIZE
        chunk["window_end"] = chunk["window_start"] + WINDOW_SIZE

        # Sum read counts within each window. The resulting methylation
        # percentage is coverage-weighted rather than an unweighted mean.
        grouped = (
            chunk.groupby(["chrom", "window_start", "window_end"], observed=True)
            .agg(
                cpg_count=("start", "size"),
                meth_sum=("meth", "sum"),
                unmeth_sum=("unmeth", "sum"),
                coverage_sum=("coverage", "sum"),
            )
            .reset_index()
        )
        summaries.append(grouped)

    if not summaries:
        raise RuntimeError(
            f"No usable CpGs found in {filepath.name}. "
            "Check the file format and coverage threshold."
        )

    result = pd.concat(summaries, ignore_index=True)

    # A genomic window may occur in more than one chunk. Merge partial
    # summaries so every window has one row per sample.
    result = (
        result.groupby(["chrom", "window_start", "window_end"], as_index=False)
        .agg(
            cpg_count=("cpg_count", "sum"),
            meth_sum=("meth_sum", "sum"),
            unmeth_sum=("unmeth_sum", "sum"),
            coverage_sum=("coverage_sum", "sum"),
        )
    )

    result["methylation"] = (
        100.0 * result["meth_sum"]
        / result["coverage_sum"].replace(0, np.nan)
    )

    result = result[result["cpg_count"] >= MIN_CPGS_PER_WINDOW].copy()
    result = result[
        ["chrom", "window_start", "window_end", "cpg_count",
         "coverage_sum", "methylation"]
    ]
    result = result.rename(
        columns={
            "cpg_count": f"{sample_name}_CpGs",
            "coverage_sum": f"{sample_name}_coverage",
            "methylation": sample_name,
        }
    )

    print(
        f"    Parsed {rows_seen:,} rows; retained {rows_kept:,} covered CpGs; "
        f"{len(result):,} windows passed the CpG-count filter.",
        flush=True,
    )
    return result


# ---------------------------------------------------------------------
# Assay-level analysis
# ---------------------------------------------------------------------

def process_assay(assay_name: str, sample_map: dict) -> None:
    print(f"\nProcessing {assay_name}...", flush=True)

    sample_tables = {}
    for sample_id, condition in sample_map.items():
        filepath = BASE / f"{sample_id}_CpG.bedGraph"
        sample_tables[condition] = summarize_one_file(filepath, condition)

    # Inner join ensures that a window is compared only when it passed
    # filtering in all four samples for this assay.
    merged = None
    for condition in ["KO_CD", "KO_HFD", "WT_CD", "WT_HFD"]:
        table = sample_tables[condition]
        if merged is None:
            merged = table
        else:
            merged = merged.merge(
                table,
                on=["chrom", "window_start", "window_end"],
                how="inner",
                validate="one_to_one",
            )

    if merged is None or merged.empty:
        raise RuntimeError(
            f"No shared windows passed filtering for {assay_name}. "
            "Check coordinate compatibility and filtering thresholds."
        )

    # Differences are percentage points: condition A minus condition B.
    merged["HFD_vs_CD_KO"] = merged["KO_HFD"] - merged["KO_CD"]
    merged["HFD_vs_CD_WT"] = merged["WT_HFD"] - merged["WT_CD"]
    merged["KO_vs_WT_CD"] = merged["KO_CD"] - merged["WT_CD"]
    merged["KO_vs_WT_HFD"] = merged["KO_HFD"] - merged["WT_HFD"]

    # Difference-in-differences: (KO HFD - KO CD) - (WT HFD - WT CD).
    merged["Genotype_by_diet_interaction"] = (
        merged["HFD_vs_CD_KO"] - merged["HFD_vs_CD_WT"]
    )

    direct_contrast_columns = [
    "HFD_vs_CD_KO",
    "HFD_vs_CD_WT",
    "KO_vs_WT_CD",
    "KO_vs_WT_HFD",
]

all_contrast_columns = direct_contrast_columns + [
    "Genotype_by_diet_interaction"
]

    # Add simple direction labels for each contrast.
    for col in contrast_columns:
        merged[f"{col}_direction"] = np.select(
            [merged[col] >= MIN_ABS_DELTA, merged[col] <= -MIN_ABS_DELTA],
            ["higher_in_A", "lower_in_A"],
            default="small_difference",
        )

    # Candidate if at least one planned contrast reaches the chosen
    # absolute-difference threshold.
    merged["max_abs_difference"] = merged[contrast_columns].abs().max(axis=1)
    candidate = merged[merged["max_abs_difference"] >= MIN_ABS_DELTA].copy()

    # Include a compact, convenient coordinate column for plotting.
    merged["region"] = (
        merged["chrom"].astype(str)
        + ":"
        + merged["window_start"].astype(str)
        + "-"
        + merged["window_end"].astype(str)
    )
    candidate["region"] = (
        candidate["chrom"].astype(str)
        + ":"
        + candidate["window_start"].astype(str)
        + "-"
        + candidate["window_end"].astype(str)
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    all_path = OUTDIR / f"{assay_name}_all_1kb_windows.tsv"
    candidate_path = OUTDIR / f"{assay_name}_candidate_DMRs.tsv"

    merged.to_csv(all_path, sep="\t", index=False)
    candidate.to_csv(candidate_path, sep="\t", index=False)

    print(f"  Shared windows: {len(merged):,}", flush=True)
    print(f"  Candidate windows: {len(candidate):,}", flush=True)
    print(f"  Wrote: {candidate_path}", flush=True)
    print(f"  Wrote: {all_path}", flush=True)


def main():
    if not BASE.exists():
        sys.exit(f"Input directory does not exist: {BASE}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    for assay_name, sample_map in ASSAYS.items():
        process_assay(assay_name, sample_map)

    print("\nFinished.", flush=True)
    print(f"Output directory: {OUTDIR}", flush=True)
    print(
        "\nNote: these are exploratory candidate regions, not statistically "
        "validated DMRs. There is one biological sample per condition.",
        flush=True,
    )


if __name__ == "__main__":
    main()
