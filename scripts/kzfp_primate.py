import pandas as pd 
import re 

# Primate specific tsv file 
primate_path = "/home/abportillo/github_repo/long-read/docs/KZFP_evolutional_primate_specific_R2.tsv"

# All KZFPs file 
kzfp_all_path = "/home/abportillo/github_repo/long-read/docs/kzfp_trono_list_fixed.csv" 

# Load data 
#Primate file 
primate_df = pd.read_csv(primate_path, sep="\t", dtype=str)
# KZFP file
all_df = pd.read_csv(kzfp_all_path, dtype=str)

# Standardize column name 
if "gene symbol" in all_df.columns:
    all_df = all_df.rename(columns={"gene symbol": "gene_symbol"})
elif "gene symbol" not in all_df.columns:
    raise ValueError("Could not find 'gene symbol or 'gene_symbol' column in kzfp_trono_list_fixed.csv")

# Clean whitespace 
primate_df["assigned_gene"] = primate_df["assigned_gene"].str.strip()
all_df["gene_symbol"] = all_df["gene_symbol"].str.strip()

#Keeping only protein-coding kzfps
primate_pc = primate_df[primate_df["classification"] == "protein_coding"].copy()

#Verify against full list
primate_set = set(primate_pc["assigned_gene"].dropna())
all_set = set(all_df["gene_symbol"].dropna())
matches = sorted(primate_set & all_set)
missing_in_all = sorted(primate_set - all_set)

# output
pd.DataFrame({"symbol": matches}).to_csv("/home/abportillo/github_repo/long-read/docs/primate_kzfp_protein_coding_matches_exact.txt", index=False, header=False)
pd.DataFrame({'symbol': missing_in_all}).to_csv("/home/abportillo/github_repo/long-read/docs/primate_kzfp_protein_coding_missing_in_full_exact.txt", index=False, header=False)