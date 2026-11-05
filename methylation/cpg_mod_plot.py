import matplotlib
matplotlib.use('Agg')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define column names
columns = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "Nvalid_cov",
    "fraction_modified", "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "Nnocall"
]

# Read the bedMethyl file
df = pd.read_csv("cpg_mod_sample2.bedMethyl", sep="\t", header=None, names=columns)

# Filter to only 5hmC and 5mC rows
df = df[df["mod_code"].isin(["h", "m"])]

# Map to modification type names
df["mod_type"] = df["mod_code"].map({"h": "5hmC", "m": "5mC"})

# Plot
plt.figure(figsize=(6, 4))
sns.boxplot(data=df, x="mod_type", y="fraction_modified", palette="Set2", showfliers=False)
sns.swarmplot(data=df.sample(n=100000, random_state=42), x="mod_type", y="fraction_modified",
              color="black", alpha=0.5, size=2)

# Style
plt.title("5hmC vs 5mC at CpG Sites")
plt.ylabel("Fraction Modified")
plt.xlabel("Cytosine Modification Type")
plt.ylim(0, 110)
plt.grid(True, axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()

# Save plot
plt.savefig("cpg_5hmC_5mC_swarm.png", dpi=300)
