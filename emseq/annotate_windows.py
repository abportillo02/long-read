#!/usr/bin/env python3
"""Annotate 10 kb windows with overlapping genes and write IGV-ready files."""
import argparse, gzip, os, re
from collections import defaultdict
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True)
ap.add_argument("--qc-dir", required=True)
ap.add_argument("--libs", default="../mapping/libraries.tsv")
ap.add_argument("--window", type=int, default=10000, help="must match qc_and_compare.py")
ap.add_argument("--max-names", type=int, default=10)
a = ap.parse_args()
W = a.window
CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
ORDER = {c: i for i, c in enumerate(CHROMS)}

def chrom_name(c):
    return "chrM" if c == "MT" else (c if c.startswith("chr") else "chr" + c)

# ---- genes -> windows ------------------------------------------------------
win_genes, n_genes = defaultdict(list), 0
op = gzip.open if a.gtf.endswith(".gz") else open
with op(a.gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        c, s, e = chrom_name(f[0]), int(f[3]) - 1, int(f[4])
        m = re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
        for w in range(s // W, (e - 1) // W + 1):
            win_genes[(c, w)].append(m.group(1))
        n_genes += 1
print(f"{n_genes:,} genes read from {a.gtf}", flush=True)

def genes_for(chrom, w):
    g = win_genes.get((chrom, w))
    if not g:
        return 0, "."
    g = sorted(set(g))
    txt = ",".join(g[:a.max_names]) + (f",+{len(g) - a.max_names}more" if len(g) > a.max_names else "")
    return len(g), txt

os.makedirs(f"{a.qc_dir}/igv", exist_ok=True)

# ---- top windows -----------------------------------------------------------
top = pd.read_csv(f"{a.qc_dir}/top_windows.tsv", sep="\t")
res = [genes_for(c, s // W) for c, s in zip(top.chrom, top.start)]
top["n_genes"], top["genes"] = [r[0] for r in res], [r[1] for r in res]
top.to_csv(f"{a.qc_dir}/top_windows.genes.tsv", sep="\t", index=False)
for assay, g in top.groupby("assay"):
    with open(f"{a.qc_dir}/igv/top_windows.{assay}.bed", "w") as fh:
        for r in g.itertuples():
            fh.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.contrast}:{r.delta:+.2f}:{r.genes}\n")
print(f"top windows: {len(top)} annotated, {(top.n_genes > 0).mean():.0%} overlap a gene")

# ---- all windows, per-library tracks and comparison tracks -----------------
meta = pd.read_csv(a.libs, sep="\t", header=None, names=["lib", "assay", "group", "ref"])

def write_bedgraph(path, name, df, col):
    with open(path, "w") as fh:
        fh.write(f'track type=bedGraph name="{name}"\n')
        df[["chrom", "start", "end", col]].to_csv(fh, sep="\t", header=False, index=False, float_format="%.3f")

for assay, gm in meta.groupby("assay"):
    info = dict(zip(gm.lib, gm.group))
    wv = pd.read_csv(f"{a.qc_dir}/window_values_{assay}.tsv.gz", sep="\t")
    wv["start"] = wv.w.astype("int64") * W
    wv["end"] = wv.start + W
    wv["o"] = wv.chrom.map(ORDER)
    wv = wv.sort_values(["o", "start"]).drop(columns="o").reset_index(drop=True)
    res = [genes_for(c, w) for c, w in zip(wv.chrom, wv.w)]
    wv["n_genes"], wv["genes"] = [r[0] for r in res], [r[1] for r in res]
    pcols = [f"{l}_pct" for l in gm.lib]
    wv[["chrom", "start", "end"] + pcols + ["n_genes", "genes"]].to_csv(
        f"{a.qc_dir}/windows_by_gene_{assay}.tsv.gz", sep="\t", index=False, float_format="%.3f")

    for l in gm.lib:
        write_bedgraph(f"{a.qc_dir}/igv/{assay}.{l}.{info[l]}.10kb_pct.bedGraph",
                       f"{assay} {l} {info[l]} %", wv, f"{l}_pct")

    tissue = [l for l in gm.lib if len(info[l].split("-")) == 3]
    def find(geno, diet):
        for l in tissue:
            p = info[l].split("-")
            if p[1] == geno and p[2] == diet:
                return l
    for name, A, B in [("KO_vs_WT_CD", find("KO", "CD"), find("WT", "CD")),
                       ("KO_vs_WT_HFD", find("KO", "HFD"), find("WT", "HFD")),
                       ("HFD_vs_CD_KO", find("KO", "HFD"), find("KO", "CD")),
                       ("HFD_vs_CD_WT", find("WT", "HFD"), find("WT", "CD"))]:
        if A and B:
            wv["delta"] = wv[f"{A}_pct"] - wv[f"{B}_pct"]
            write_bedgraph(f"{a.qc_dir}/igv/{assay}.{name}.delta.bedGraph",
                           f"{assay} {name} (A minus B, % points)", wv, "delta")
    print(f"{assay}: {len(wv):,} windows, {(wv.n_genes > 0).mean():.0%} overlap a gene")
print("done")