#!/usr/bin/env python3

import glob
import os
import pandas as pd

INPUT_DIR = "/net/nfs-irwrsrchnas01/labs/yunroseli_grp/mchosco/Abner/emseq-test/mapping/mbias"
OUTPUT_FILE = os.path.join(INPUT_DIR, "mbias_summary.tsv")

# Define a central plateau region.
# Review the plots to confirm this interval is appropriate.
PLATEAU_START = 10
PLATEAU_END = 140

FIRST_N = 6
LAST_N = 6

rows = []

for path in sorted(glob.glob(os.path.join(INPUT_DIR, "*.mbias.tsv"))):
    library = os.path.basename(path).replace(".mbias.tsv", "")

    df = pd.read_csv(path, sep="\t")

    required = {
        "Strand", "Read", "Position",
        "nMethylated", "nUnmethylated"
    }
    if not required.issubset(df.columns):
        raise ValueError(
            f"Unexpected columns in {path}: {list(df.columns)}"
        )

    df["Read"] = df["Read"].astype(int)
    df["Position"] = df["Position"].astype(int)

    total = df["nMethylated"] + df["nUnmethylated"]
    df["Methylation_pct"] = (
        100 * df["nMethylated"] / total.where(total > 0)
    )

    for (strand, read), group in df.groupby(["Strand", "Read"]):
        group = group.sort_values("Position")
        valid = group.dropna(subset=["Methylation_pct"])

        first = valid[valid["Position"] <= FIRST_N]
        last_positions = sorted(valid["Position"].unique())[-LAST_N:]
        last = valid[valid["Position"].isin(last_positions)]

        plateau = valid[
            valid["Position"].between(PLATEAU_START, PLATEAU_END)
        ]

        rows.append({
            "library": library,
            "strand": strand,
            "read": f"R{read}",
            "plateau_pct": plateau["Methylation_pct"].mean(),
            "first6_pct": first["Methylation_pct"].mean(),
            "last6_pct": last["Methylation_pct"].mean(),
            "plateau_positions": len(plateau),
            "first_positions": len(first),
            "last_positions": len(last),
            "mean_plateau_coverage": (
                plateau["nMethylated"] + plateau["nUnmethylated"]
            ).mean(),
            "max_position": (
                int(valid["Position"].max()) if not valid.empty else None
            ),
        })

summary = pd.DataFrame(rows)
summary.to_csv(OUTPUT_FILE, sep="\t", index=False, float_format="%.3f")

print(f"Wrote {len(summary)} rows to {OUTPUT_FILE}")
print(summary.to_string(index=False))
