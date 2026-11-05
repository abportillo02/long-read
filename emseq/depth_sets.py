#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--libs", default="../mapping/libraries.tsv")
ap.add_argument("--dir", required=True)
ap.add_argument("--pattern", default="{lib}.grcm39.counts.bed")  # cols: chr start end meth unmeth
ap.add_argument("--thresholds", type=int, nargs="+", default=[1, 3, 5, 10, 15, 20, 30])
ap.add_argument("--out", default="depth_sets.tsv")
a = ap.parse_args()

meta = pd.read_csv(a.libs, sep="\t", header=None, names=["lib", "assay", "group", "ref"])
depth = {}                       # lib -> {chrom: (starts, depths)}
for lib in meta.lib:
    df = pd.read_csv(f"{a.dir}/{a.pattern.format(lib=lib)}", sep="\t", header=None,
                     usecols=[0, 1, 3, 4], names=["c", "s", "m", "u"],
                     dtype={"c": "category", "s": "int32", "m": "int32", "u": "int32"})
    df["d"] = df.m + df.u
    depth[lib] = {c: (g.s.to_numpy(), g.d.to_numpy()) for c, g in df.groupby("c", observed=True)}
    print(f"{lib}: {len(df):,} CpGs, median depth {int(df.d.median())}", flush=True)

sets = {"all_10": list(meta.lib)}
for assay, g in meta.groupby("assay"):
    sets[f"{assay}_all"] = list(g.lib)
    sets[f"{assay}_no_CT26"] = list(g[g.group != "CT26"].lib)

rows = []
for name, libs in sets.items():
    chroms = set.intersection(*[set(depth[l]) for l in libs])
    for t in a.thresholds:
        n = 0
        for c in chroms:
            keep = None
            for l in libs:
                s, d = depth[l][c]
                p = s[d >= t]
                keep = p if keep is None else np.intersect1d(keep, p, assume_unique=True)
            n += len(keep)
        rows.append((name, len(libs), t, n))
        print(name, len(libs), t, f"{n:,}", flush=True)
pd.DataFrame(rows, columns=["set", "n_libs", "min_depth", "shared_cpgs"]).to_csv(a.out, sep="\t", index=False)