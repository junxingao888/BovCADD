import pandas as pd
import glob
import sys

files = glob.glob("results/burden/chr*.hom_burden.csv")
dfs = [pd.read_csv(f) for f in files]

merged = dfs[0][["sample_id"]].copy()
for col in ["hom_deleterious", "het_deleterious", "hom_benign"]:
    merged[col] = sum(df[col] for df in dfs)

merged["burden_ratio"] = merged["hom_deleterious"] / merged["hom_benign"].replace(0, pd.NA)
merged.to_csv("results/burden/genomewide_hom_burden.csv", index=False)
print("Wrote results/burden/genomewide_hom_burden.csv")
