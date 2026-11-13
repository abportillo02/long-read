#!/usr/bin/env python3
"""
Methylation heatmaps and individual gene-locus plots for the AMPK mouse colon
EM-seq and 5hmC-seq data.

Expected inputs (see the settings block below):
  mapping/methyl/grcm39_v2/YRLP_00NN.grcm39.counts.bed   (lifted to GRCm39)
  mapping/DMR_analysis/EM_seq_candidate_DMRs.tsv
  mapping/DMR_analysis/5hmC_seq_candidate_DMRs.tsv
  annotation/gencode.vM39.annotation.gtf

Outputs:
  mapping/DMR_analysis/figures/
    EM_seq_heatmap_absolute.png
    EM_seq_heatmap_row_zscore.png
    5hmC_seq_heatmap_absolute.png
    5hmC_seq_heatmap_row_zscore.png
    loci/<assay>_<gene>_methylation_locus.png

Notes
  * CpG coordinates, GTF coordinates and candidate windows must all be in
    GRCm39 coordinates. The unlifted bedGraphs in mapping/methyl/ are in
    sample-specific reference coordinates and must NOT be mixed with the GTF.
  * One biological sample per condition: plots are descriptive, not
    replicate-based statistics.
  * Each counts file is read once per assay, for all genes together.
"""

from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# PATHS AND SETTINGS (edit these)
# ---------------------------------------------------------------------

PROJECT = Path(
    "/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test"
)

# Lifted (GRCm39) per-CpG counts. To plot from the depth-filtered tables
# instead, point COUNTS_DIR at .../grcm39_v2/filtered_d5, set COUNTS_PATTERN
# to the real file name there, and set COUNT_COLUMNS to the real column order.
COUNTS_DIR = PROJECT / "mapping" / "methyl" / "grcm39_v2"
COUNTS_PATTERN = "{sample}.grcm39.counts.bed"
# Column order in the counts file: chrom, start (0-based), end, then the
# methylated and unmethylated read counts. Only the first five columns are read.
COUNT_COLUMNS = ["chrom", "start", "end", "meth", "unmeth"]

DMR_DIR = PROJECT / "mapping" / "DMR_analysis_grcm39_v2_filtered_d5"
GTF = PROJECT / "annotation" / "gencode.vM39.annotation.gtf"
FIGURE_DIR = DMR_DIR / "figures"

# Gene symbols (mouse nomenclature, as in the GTF gene_name field).
GENES_OF_INTEREST = ["Cdx2"]

CHUNK_SIZE = 2_000_000
MIN_CPG_COVERAGE = 5      # per-CpG depth threshold used for the locus dots
FLANK_BP = 20_000
TOP_N_HEATMAP = 100

# Upper end of the color scale (heatmap) and y-axis (locus plot), in percent.
# 5hmC levels are low (genome-wide means are roughly 2-3%), so a 0-100 scale
# would flatten them. Adjust after looking at the first figures.
PLOT_MAX = {"EM_seq": 100, "5hmC_seq": 20}

ASSAYS = {
    "EM_seq": {
        "KO_CD": "YRLP_0064",
        "KO_HFD": "YRLP_0069",
        "WT_CD": "YRLP_0075",
        "WT_HFD": "YRLP_0065",
    },
    "5hmC_seq": {
        "KO_CD": "YRLP_0076",
        "KO_HFD": "YRLP_0072",
        "WT_CD": "YRLP_0074",
        "WT_HFD": "YRLP_0071",
    },
}

CONDITIONS = ["KO_CD", "KO_HFD", "WT_CD", "WT_HFD"]
COLORS = {
    "KO_CD": "#3366cc",
    "KO_HFD": "#dc3912",
    "WT_CD": "#109618",
    "WT_HFD": "#990099",
}

# Set at runtime by preflight(): whether the counts files use "chr1" style
# names (True) or "1" style names (False).
COUNTS_HAS_CHR = True


# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------

def strip_chr(name):
    """Remove a leading 'chr' so chromosome names can be compared."""
    name = str(name)
    return name[3:] if name.lower().startswith("chr") else name


def counts_chrom_name(gtf_chrom):
    """Convert a GTF-style chromosome name to the style in the counts files."""
    return "chr" + strip_chr(gtf_chrom) if COUNTS_HAS_CHR else strip_chr(gtf_chrom)


def counts_path(sample_id):
    return COUNTS_DIR / COUNTS_PATTERN.format(sample=sample_id)


def detect_skiprows(path, max_lines=10):
    """Count leading header lines (e.g. 'track ...') before the first data row.

    A data row is any line whose second tab-separated field is an integer.
    """
    skip = 0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2 and fields[1].isdigit():
                break
            skip += 1
            if skip >= max_lines:
                break
    return skip


def peek_data_lines(path, n=2):
    skip = detect_skiprows(path)
    lines = []
    with open(path) as handle:
        for i, line in enumerate(handle):
            if i < skip:
                continue
            lines.append(line.rstrip("\n"))
            if len(lines) >= n:
                break
    return lines


def parse_gtf_attributes(attribute_text):
    """Parse GTF attributes such as gene_name and gene_id."""
    attrs = {}
    for item in attribute_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'(\S+)\s+"([^"]*)"', item)
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs


def load_gene_annotation(gtf_path, genes_of_interest):
    """Load gene coordinates and exon coordinates for requested symbols."""
    wanted = set(genes_of_interest)
    gene_records = {}
    exon_records = []

    print(f"Reading GTF annotation: {gtf_path}", flush=True)
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr_text = fields
            if feature not in {"gene", "exon"}:
                continue

            attrs = parse_gtf_attributes(attr_text)
            gene_name = attrs.get("gene_name")
            if gene_name not in wanted:
                continue

            # GTF is 1-based inclusive. Convert to 0-based, half-open.
            record = {
                "gene": gene_name,
                "gene_id": attrs.get("gene_id", ""),
                "chrom": chrom,
                "start": int(start) - 1,
                "end": int(end),
                "strand": strand,
            }

            if feature == "gene":
                prior = gene_records.get(gene_name)
                if prior is not None:
                    print(
                        f"NOTE: {gene_name} has more than one gene feature "
                        f"in the GTF; keeping the longer one.",
                        flush=True,
                    )
                    if (prior["end"] - prior["start"]) >= (
                        record["end"] - record["start"]
                    ):
                        continue
                gene_records[gene_name] = record
            else:
                exon_records.append(record)

    missing = sorted(wanted - set(gene_records))
    if missing:
        print(
            "WARNING: These gene symbols were not found as gene features in "
            f"the GTF: {', '.join(missing)}",
            flush=True,
        )

    return gene_records, pd.DataFrame(exon_records)


def read_cpgs_for_regions(filepath, regions):
    """Read one counts file once and return covered CpGs for every region.

    regions: dict gene -> (chrom in GTF naming, region_start, region_end)
    returns: dict gene -> DataFrame[start, methylation, coverage]
    """
    wanted = {
        gene: (counts_chrom_name(chrom), start, end)
        for gene, (chrom, start, end) in regions.items()
    }
    needed_chroms = {chrom for chrom, _, _ in wanted.values()}
    pieces = {gene: [] for gene in regions}

    reader = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        skiprows=detect_skiprows(filepath),
        names=COUNT_COLUMNS,
        usecols=range(len(COUNT_COLUMNS)),
        dtype={
            "chrom": str, "start": "int64", "end": "int64",
            "meth": "float64", "unmeth": "float64",
        },
        comment="#",
        chunksize=CHUNK_SIZE,
        on_bad_lines="skip",
    )

    for chunk in reader:
        chunk = chunk[chunk["chrom"].isin(needed_chroms)]
        if chunk.empty:
            continue
        coverage = chunk["meth"] + chunk["unmeth"]
        keep = (coverage >= MIN_CPG_COVERAGE) & (chunk["meth"] >= 0) & (
            chunk["unmeth"] >= 0
        )
        chunk = chunk[keep]
        coverage = coverage[keep]
        if chunk.empty:
            continue

        for gene, (chrom, start, end) in wanted.items():
            mask = (
                (chunk["chrom"] == chrom)
                & (chunk["start"] < end)
                & (chunk["end"] > start)
            )
            if not mask.any():
                continue
            pieces[gene].append(
                pd.DataFrame({
                    "start": chunk.loc[mask, "start"],
                    "methylation": 100.0 * chunk.loc[mask, "meth"] / coverage[mask],
                    "coverage": coverage[mask],
                })
            )

    result = {}
    for gene, parts in pieces.items():
        if parts:
            result[gene] = pd.concat(parts, ignore_index=True).sort_values("start")
        else:
            result[gene] = pd.DataFrame(columns=["start", "methylation", "coverage"])
    return result


# ---------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------

def preflight():
    """Fail early on missing files and print what to eyeball."""
    global COUNTS_HAS_CHR

    problems = []
    if not GTF.exists():
        problems.append(f"GTF not found: {GTF}")
    for assay, samples in ASSAYS.items():
        for condition, sample_id in samples.items():
            path = counts_path(sample_id)
            if not path.exists():
                problems.append(f"Counts file not found ({assay} {condition}): {path}")
    if problems:
        sys.exit("\n".join(["PREFLIGHT FAILED:"] + problems))

    first = counts_path(next(iter(ASSAYS["EM_seq"].values())))
    lines = peek_data_lines(first, n=2)
    print(f"First data lines of {first.name}:", flush=True)
    for line in lines:
        print(f"  {line}", flush=True)
    COUNTS_HAS_CHR = lines[0].split("\t")[0].lower().startswith("chr")
    print(
        f"Counts files use {'chr-prefixed' if COUNTS_HAS_CHR else 'unprefixed'} "
        "chromosome names; GTF names are converted to match.",
        flush=True,
    )

    for assay in ASSAYS:
        dmr_file = DMR_DIR / f"{assay}_candidate_DMRs.tsv"
        if dmr_file.exists():
            head = pd.read_csv(dmr_file, sep="\t", nrows=2)
            print(f"{dmr_file.name} columns: {list(head.columns)}", flush=True)
            print(head.head(2).to_string(index=False), flush=True)
        else:
            print(f"WARNING: candidate table not found: {dmr_file}", flush=True)
    print("", flush=True)


# ---------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------

def make_heatmaps(assay_name):
    dmr_file = DMR_DIR / f"{assay_name}_candidate_DMRs.tsv"
    if not dmr_file.exists():
        print(f"Skipping heatmap; candidate table not found: {dmr_file}")
        return

    df = pd.read_csv(dmr_file, sep="\t")
    required = ["chrom", "window_start", "window_end"] + CONDITIONS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{dmr_file.name} is missing required columns: {missing}"
        )

    df = df.dropna(subset=CONDITIONS).copy()
    if df.empty:
        print(f"No complete candidate rows in {dmr_file}")
        return

    # If the table stores fractions (0-1) instead of percent, convert.
    scale = 1.0
    if df[CONDITIONS].to_numpy(dtype=float).max() <= 1.0:
        scale = 100.0
        print(
            f"WARNING: {dmr_file.name} values look like fractions (max <= 1); "
            "multiplying by 100 for plotting.",
            flush=True,
        )
        df[CONDITIONS] = df[CONDITIONS] * scale

    # Rank windows by the largest difference across the planned contrasts.
    if "max_abs_difference" in df.columns:
        df = df.sort_values("max_abs_difference", ascending=False)
    else:
        values = df[CONDITIONS]
        df["_range"] = values.max(axis=1) - values.min(axis=1)
        df = df.sort_values("_range", ascending=False)

    plot_df = df.head(TOP_N_HEATMAP).copy()
    labels = (
        plot_df["chrom"].astype(str)
        + ":"
        + plot_df["window_start"].astype(int).astype(str)
        + "-"
        + plot_df["window_end"].astype(int).astype(str)
    )
    matrix = plot_df[CONDITIONS].to_numpy(dtype=float)
    vmax = PLOT_MAX[assay_name]

    # Absolute methylation heatmap.
    fig_height = max(6, 0.22 * len(plot_df))
    fig, ax = plt.subplots(figsize=(8, fig_height))
    im = ax.imshow(
        matrix, aspect="auto", interpolation="nearest",
        cmap="viridis", vmin=0, vmax=vmax
    )
    ax.set_xticks(range(len(CONDITIONS)))
    ax.set_xticklabels(CONDITIONS, rotation=35, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title(f"{assay_name}: top {len(plot_df)} candidate windows")
    ax.set_ylabel("Genomic window")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"Methylation (%), color scale 0-{vmax}")
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{assay_name}_heatmap_absolute.png",
        dpi=250, bbox_inches="tight"
    )
    plt.close(fig)

    # Row z-score heatmap emphasizes relative condition patterns per region.
    row_mean = matrix.mean(axis=1, keepdims=True)
    row_sd = matrix.std(axis=1, keepdims=True)
    z = np.divide(
        matrix - row_mean,
        row_sd,
        out=np.zeros_like(matrix),
        where=row_sd > 0,
    )
    fig, ax = plt.subplots(figsize=(8, fig_height))
    im = ax.imshow(
        z, aspect="auto", interpolation="nearest",
        cmap="RdBu_r", vmin=-2, vmax=2
    )
    ax.set_xticks(range(len(CONDITIONS)))
    ax.set_xticklabels(CONDITIONS, rotation=35, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title(f"{assay_name}: row-scaled methylation patterns")
    ax.set_ylabel("Genomic window")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Row Z-score")
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{assay_name}_heatmap_row_zscore.png",
        dpi=250, bbox_inches="tight"
    )
    plt.close(fig)

    print(
        f"Saved heatmaps for {assay_name} ({len(plot_df)} regions).",
        flush=True,
    )


# ---------------------------------------------------------------------
# Individual gene-region plots
# ---------------------------------------------------------------------

def make_gene_locus_plot(assay_name, gene, gene_record, exons, cpg_tables):
    chrom = gene_record["chrom"]
    gene_start = gene_record["start"]
    gene_end = gene_record["end"]
    strand = gene_record["strand"]

    region_start = max(0, gene_start - FLANK_BP)
    region_end = gene_end + FLANK_BP
    ymax = PLOT_MAX[assay_name]

    # Candidate windows overlapping this plotted locus (chromosome names are
    # compared without any 'chr' prefix).
    dmr_file = DMR_DIR / f"{assay_name}_candidate_DMRs.tsv"
    candidate = pd.DataFrame()
    if dmr_file.exists():
        dmr = pd.read_csv(dmr_file, sep="\t")
        candidate = dmr[
            (dmr["chrom"].map(strip_chr) == strip_chr(chrom))
            & (dmr["window_start"] < region_end)
            & (dmr["window_end"] > region_start)
        ].copy()

    fig, axes = plt.subplots(
        nrows=5, ncols=1, figsize=(13, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.8]},
    )

    for ax, condition in zip(axes[:4], CONDITIONS):
        sample_id = ASSAYS[assay_name][condition]
        cpgs = cpg_tables[condition][gene]

        if cpgs.empty:
            print(
                f"WARNING: no CpGs (depth >= {MIN_CPG_COVERAGE}) for "
                f"{sample_id} ({assay_name} {condition}) at {gene}. Check "
                "chromosome naming and the counts file.",
                flush=True,
            )
        else:
            ax.scatter(
                cpgs["start"] / 1e6,
                cpgs["methylation"].clip(upper=ymax),
                s=9,
                alpha=0.65,
                color=COLORS[condition],
                edgecolors="none",
            )

        # Highlight candidate windows behind the CpG points.
        for _, row in candidate.iterrows():
            ax.axvspan(
                row["window_start"] / 1e6,
                row["window_end"] / 1e6,
                color="gold",
                alpha=0.16,
                linewidth=0,
            )

        ax.set_ylim(-0.03 * ymax, 1.03 * ymax)
        ax.set_ylabel(f"{condition}\n{sample_id}\n(%)", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Gene model: gene line plus exon boxes (exons of this gene_id only).
    gene_ax = axes[4]
    gene_ax.hlines(0.5, gene_start / 1e6, gene_end / 1e6,
                   color="black", linewidth=1.5)

    if exons.empty:
        gene_exons = pd.DataFrame()
    else:
        gene_exons = exons[
            (exons["gene_id"] == gene_record["gene_id"])
            & (exons["chrom"] == chrom)
        ]

    for _, exon in gene_exons.iterrows():
        gene_ax.add_patch(
            plt.Rectangle(
                (exon["start"] / 1e6, 0.32),
                (exon["end"] - exon["start"]) / 1e6,
                0.36,
                color="black",
            )
        )

    direction = "→" if strand == "+" else "←"
    gene_ax.text(
        (gene_start + gene_end) / 2 / 1e6,
        0.82,
        f"{gene} ({direction})",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    gene_ax.set_xlim(region_start / 1e6, region_end / 1e6)
    gene_ax.set_ylim(0, 1.25)
    gene_ax.set_yticks([])
    gene_ax.set_xlabel(f"{chrom} genomic position (Mb, GRCm39)")
    gene_ax.spines["left"].set_visible(False)
    gene_ax.spines["right"].set_visible(False)
    gene_ax.spines["top"].set_visible(False)

    axes[0].set_title(
        f"{assay_name}: {gene} locus | "
        f"{chrom}:{region_start:,}-{region_end:,} (GRCm39)"
    )
    fig.suptitle(
        f"Dots = individual CpGs (depth ≥ {MIN_CPG_COVERAGE}, y-axis 0-{ymax}%); "
        "gold shading = candidate windows",
        y=0.995,
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    locus_dir = FIGURE_DIR / "loci"
    locus_dir.mkdir(parents=True, exist_ok=True)
    out = locus_dir / f"{assay_name}_{gene}_methylation_locus.png"
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}", flush=True)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    preflight()

    gene_records, exons = load_gene_annotation(GTF, GENES_OF_INTEREST)
    if not gene_records:
        sys.exit(
            "None of the requested genes were found. Check the gene symbols "
            "in GENES_OF_INTEREST and the GTF annotation."
        )

    regions = {
        gene: (
            rec["chrom"],
            max(0, rec["start"] - FLANK_BP),
            rec["end"] + FLANK_BP,
        )
        for gene, rec in gene_records.items()
    }

    for assay_name, samples in ASSAYS.items():
        print(f"\nMaking {assay_name} heatmaps...", flush=True)
        make_heatmaps(assay_name)

        # Read each counts file once for all genes of interest.
        cpg_tables = {}
        for condition in CONDITIONS:
            path = counts_path(samples[condition])
            print(f"Reading {path.name} ({assay_name} {condition})...", flush=True)
            cpg_tables[condition] = read_cpgs_for_regions(path, regions)

        for gene, record in gene_records.items():
            print(f"Plotting {assay_name} locus: {gene}", flush=True)
            make_gene_locus_plot(assay_name, gene, record, exons, cpg_tables)

    print("\nFinished creating figures.", flush=True)
    print(f"Figures saved under: {FIGURE_DIR}", flush=True)


if __name__ == "__main__":
    main()