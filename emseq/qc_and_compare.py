#!/usr/bin/env python3
"""QC figures and simple 10 kb window comparison from the filtered CpG tables."""
import argparse, os, shutil, datetime
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--libs", default="../mapping/libraries.tsv")
ap.add_argument("--dir", required=True, help="folder with <lib>.filtered.tsv.gz and filter_summary.tsv")
ap.add_argument("--outdir", required=True)
ap.add_argument("--window", type=int, default=10000)
ap.add_argument("--min-cpgs", type=int, default=10)
ap.add_argument("--min-calls", type=int, default=200)
ap.add_argument("--top", type=int, default=200)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
ORDER = {c: i for i, c in enumerate(CHROMS)}
meta = pd.read_csv(a.libs, sep="\t", header=None, names=["lib", "assay", "group", "ref"])
summ = pd.read_csv(f"{a.dir}/filter_summary.tsv", sep="\t").set_index("library")
info_all = dict(zip(meta.lib, meta.group))

# ---- aggregate each library into windows ----------------------------------
win, hist = {}, {}
for l in meta.lib:
    df = pd.read_csv(f"{a.dir}/{l}.filtered.tsv.gz", sep="\t",
                     usecols=["chrom", "start", "n_meth", "n_unmeth"],
                     dtype={"chrom": "category", "start": "int32", "n_meth": "int32", "n_unmeth": "int32"})
    hist[l] = np.bincount(np.minimum((df.n_meth + df.n_unmeth).to_numpy(), 120), minlength=121)
    df["w"] = (df.start // a.window).astype("int32")
    g = df.groupby(["chrom", "w"], observed=True).agg(
        m=("n_meth", "sum"), u=("n_unmeth", "sum"), n=("start", "size")).astype("int64")
    g = g[(g.n >= a.min_cpgs) & ((g.m + g.u) >= a.min_calls)].copy()
    g["pct"] = 100 * g.m / (g.m + g.u)
    win[l] = g
    print(f"{l}: {len(g):,} windows", flush=True)

def label(l): return f"{l[-4:]}\n{info_all[l]}"

# ---- figure 1: global means and depth -------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for i, (assay, gm) in enumerate(meta.groupby("assay")):
    ls = list(gm.sort_values("group").lib)
    ax[i].bar(range(len(ls)), [summ.loc[l, "mean_meth_pct"] for l in ls], color="steelblue")
    ax[i].set_xticks(range(len(ls))); ax[i].set_xticklabels([label(l) for l in ls], fontsize=7)
    ax[i].set_ylabel("mean methylation (%)"); ax[i].set_title(f"{assay}: global mean")
for l in meta.lib:
    ax[2].plot(np.arange(121), hist[l] / hist[l].sum(), lw=1, label=l[-4:])
ax[2].set_xlabel("depth per CpG (capped at 120)"); ax[2].set_ylabel("fraction of CpGs")
ax[2].set_title("Depth after filtering"); ax[2].legend(fontsize=6, ncol=2)
plt.tight_layout(); plt.savefig(f"{a.outdir}/fig1_global_and_depth.png", dpi=150); plt.close()

# ---- per assay: heatmap, contrasts, scatter -------------------------------
summary_rows, top_all = [], []
for assay, gm in meta.groupby("assay"):
    libs = list(gm.lib)
    info = dict(zip(gm.lib, gm.group))
    tissue = [l for l in libs if len(info[l].split("-")) == 3]

    # correlation heatmap, all libraries of the assay
    pct = pd.concat({l: win[l]["pct"] for l in libs}, axis=1, join="inner")
    corr = pct.corr().to_numpy()
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    im = ax.imshow(corr, vmin=min(0.5, corr.min()), vmax=1, cmap="viridis")
    ax.set_xticks(range(len(libs))); ax.set_yticks(range(len(libs)))
    ax.set_xticklabels([label(l) for l in libs], fontsize=6); ax.set_yticklabels([label(l) for l in libs], fontsize=6)
    for i in range(len(libs)):
        for j in range(len(libs)):
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", color="white", fontsize=6)
    plt.colorbar(im, ax=ax, label="Pearson r (10 kb windows)"); ax.set_title(f"{assay}: library correlation")
    plt.tight_layout(); plt.savefig(f"{a.outdir}/fig2_corr_{assay}.png", dpi=150); plt.close()
    out = pct.copy(); out.columns = [f"{c}_pct" for c in libs]
    out.reset_index().to_csv(f"{a.outdir}/window_values_{assay}.tsv.gz", sep="\t", index=False)

    # contrasts between tissue libraries
    def find(geno, diet):
        for l in tissue:
            p = info[l].split("-")
            if p[1] == geno and p[2] == diet:
                return l
    contrasts = [("KO_vs_WT_CD", find("KO", "CD"), find("WT", "CD")),
                 ("KO_vs_WT_HFD", find("KO", "HFD"), find("WT", "HFD")),
                 ("HFD_vs_CD_KO", find("KO", "HFD"), find("KO", "CD")),
                 ("HFD_vs_CD_WT", find("WT", "HFD"), find("WT", "CD"))]
    idx = win[tissue[0]].index
    for l in tissue[1:]:
        idx = idx.intersection(win[l].index)
    W = {l: win[l].loc[idx] for l in tissue}
    chrom = idx.get_level_values(0).astype(str)
    start = idx.get_level_values(1).to_numpy().astype("int64") * a.window

    fig, axes = plt.subplots(2, 2, figsize=(10, 9)); axes = axes.ravel()
    for k, (name, A, B) in enumerate(contrasts):
        if A is None or B is None:
            continue
        mA, uA, mB, uB = W[A].m.to_numpy(), W[A].u.to_numpy(), W[B].m.to_numpy(), W[B].u.to_numpy()
        nA, nB = mA + uA, mB + uB
        pA, pB = mA / nA, mB / nB
        p = (mA + mB) / (nA + nB)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (pA - pB) / np.sqrt(p * (1 - p) * (1 / nA + 1 / nB))
        z = np.where(np.isfinite(z), z, 0.0)
        gA, gB = 100 * mA.sum() / nA.sum(), 100 * mB.sum() / nB.sum()
        r = float(np.corrcoef(pA, pB)[0, 1])
        summary_rows.append(dict(assay=assay, contrast=name, A=f"{A} {info[A]}", B=f"{B} {info[B]}",
                                 mean_A_pct=round(gA, 3), mean_B_pct=round(gB, 3),
                                 diff_A_minus_B=round(gA - gB, 3), window_corr=round(r, 4), n_windows=len(idx)))
        top = pd.DataFrame({"assay": assay, "contrast": name, "chrom": chrom, "start": start, "end": start + a.window,
                            "pct_A": 100 * pA, "pct_B": 100 * pB, "delta": 100 * (pA - pB), "z": z,
                            "calls_A": nA, "calls_B": nB})
        top_all.append(top.reindex(top.z.abs().sort_values(ascending=False).index).head(a.top))
        hi = float(np.percentile(np.r_[100 * pA, 100 * pB], 99.9)) * 1.05
        axes[k].hexbin(100 * pB, 100 * pA, gridsize=70, bins="log", mincnt=1, cmap="viridis", extent=(0, hi, 0, hi))
        axes[k].plot([0, hi], [0, hi], "r--", lw=0.8)
        axes[k].set_xlabel(f"{B[-4:]} {info[B]} (%)"); axes[k].set_ylabel(f"{A[-4:]} {info[A]} (%)")
        axes[k].set_title(f"{name}\nr={r:.3f}; means {gA:.2f}% vs {gB:.2f}% (diff {gA - gB:+.2f})", fontsize=9)
    plt.suptitle(f"{assay}: 10 kb window methylation"); plt.tight_layout()
    plt.savefig(f"{a.outdir}/fig3_compare_{assay}.png", dpi=150); plt.close()

sm = pd.DataFrame(summary_rows)
sm.to_csv(f"{a.outdir}/group_comparison_summary.tsv", sep="\t", index=False)
tp = pd.concat(top_all)
tp["o"] = tp.chrom.map(ORDER)
tp.sort_values(["assay", "contrast", "o", "start"]).drop(columns="o").to_csv(
    f"{a.outdir}/top_windows.tsv", sep="\t", index=False)
with open(f"{a.outdir}/qc_params.txt", "w") as fh:
    fh.write(f"date: {datetime.datetime.now().isoformat()}\npandas {pd.__version__}, numpy {np.__version__}, "
             f"matplotlib {matplotlib.__version__}\nargs: {vars(a)}\n"
             "z is a two-proportion ranking aid on pooled reads, not a p-value (one library per condition).\n")
shutil.copy(os.path.abspath(__file__), a.outdir)
print(sm.to_string(index=False))