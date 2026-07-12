import pandas as pd
from scipy.stats import mannwhitneyu
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

burden = pd.read_csv("results/burden/genomewide_hom_burden.csv")
meta = pd.read_csv("config/sample_breed_map.csv")  # sample_id, breed
df = burden.merge(meta, on="sample_id")

# Pairwise comparisons of burden_ratio between breeds
results = []
for b1, b2 in itertools.combinations(df.breed.unique(), 2):
    g1 = df[df.breed == b1].burden_ratio.dropna()
    g2 = df[df.breed == b2].burden_ratio.dropna()
    stat, p = mannwhitneyu(g1, g2)
    results.append({"breed1": b1, "breed2": b2, "p_value": p})

pd.DataFrame(results).to_csv("results/burden_summary.csv", index=False)

# Plot
plt.figure(figsize=(8, 5))
sns.boxplot(data=df, x="breed", y="burden_ratio")
plt.ylabel("Hom(BovCADD>20) / Hom(BovCADD<5)")
plt.title("Genetic burden by breed (BovCADD homozygosity ratio)")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("results/plots/burden_by_breed.png", dpi=300)
