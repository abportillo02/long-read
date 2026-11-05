#!/usr/bin/env python3
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="../mapping/mbias")
ap.add_argument("--libs", nargs="+", required=True)
ap.add_argument("--scan", type=int, default=50)
ap.add_argument("--tol-abs", type=float, default=0.4)
ap.add_argument("--tol-rel", type=float, default=0.03)
ap.add_argument("--out", default="edge_profile.tsv")
a = ap.parse_args()

profiles, trims = [], []
for lib in a.libs:
    df = pd.read_csv(f"{a.dir}/{lib}.mbias.tsv", sep="\t")
    df["cov"] = df.nMethylated + df.nUnmethylated
    df["pct"] = 100 * df.nMethylated / df["cov"]
    for (strand, read), g in df.groupby(["Strand", "Read"]):
        g = g.sort_values("Position")
        n = int(g.Position.max())
        core = g[(g.Position > a.scan) & (g.Position <= n - a.scan)]
        plateau = core.pct.median()
        tol = max(a.tol_abs, a.tol_rel * plateau)
        g = g.assign(library=lib, plateau=plateau, delta=g.pct - plateau)
        profiles.append(g)
        bad = g[g.delta.abs() > tol]
        bl = bad[bad.Position <= a.scan].Position
        br = bad[bad.Position > n - a.scan].Position
        left = int(bl.max()) if len(bl) else 0
        right = int(n - br.min() + 1) if len(br) else 0
        trims.append(dict(library=lib, strand=strand, read=read, n=n,
                          plateau=round(plateau, 2), tol=round(tol, 2),
                          left_trim=left, right_trim=right,
                          capped=(left >= a.scan or right >= a.scan)))

pd.concat(profiles).to_csv(a.out, sep="\t", index=False)
t = pd.DataFrame(trims)
print(t.to_string(index=False))

print("\nMost conservative trim across libraries (alignment coordinates):")
m = t.groupby(["strand", "read"]).agg(n=("n", "max"), left=("left_trim", "max"),
                                       right=("right_trim", "max")).reset_index()
print(m.to_string(index=False))

def bounds(strand):
    out = []
    for read in (1, 2):
        r = m[(m.strand == strand) & (m.read == read)].iloc[0]
        out.append(r.left + 1 if r.left > 0 else 0)
        out.append(r.n - r.right if r.right > 0 else 0)
    return ",".join(str(int(x)) for x in out)

print(f"\nSuggested common window: --OT {bounds('OT')} --OB {bounds('OB')}")

# python3 edge_trim.py --dir ../mapping/mbias \
#   --libs YRLP_0071 YRLP_0072 YRLP_0073 YRLP_0074 YRLP_0076 \
#   --out edge_profile_5hmC.tsv

# python3 edge_trim.py --dir ../mapping/mbias \
#   --libs YRLP_0064 YRLP_0065 YRLP_0069 YRLP_0070 YRLP_0075 \
#   --out edge_profile_EMseq.tsv

# Table to look at residuals

# python3 - <<'EOF'
# import pandas as pd
# d = pd.read_csv("edge_profile_EMseq.tsv", sep="\t")
# x = d[(d.Strand=="OB") & (d.Read==2) & (d.Position<=25)]
# print(x.pivot(index="Position", columns="library", values="delta").round(2).to_string())
# h = pd.read_csv("edge_profile_5hmC.tsv", sep="\t")
# y = h[(h.Strand=="OB") & (h.Read==2) & (h.Position>=118)]
# print(y.pivot(index="Position", columns="library", values="delta").round(2).to_string())
# EOF
