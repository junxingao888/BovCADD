# Cattle Genetic Burden Pipeline (BovCADD-based)

A Linux/Python pipeline to estimate **genetic burden** in cattle populations
using **BovCADD** deleteriousness scores, developed for the 1000 Bull Genomes
(Run9) OPTIBOV ~4000 TAUIND dataset. Metrics are computed as **percentages of genotyped
sites** to control for breed admixture (taurine × indicine ancestry), since
raw variant counts are confounded by ancestry composition and site
ascertainment differences across breeds.

---

## 1. Background

Individuals from admixed or genetically diverse cattle populations differ
substantially in the *number* of variable/callable sites relative to a
reference panel. Comparing raw homozygous/heterozygous counts across breeds
is therefore misleading. To correct for this, burden is expressed as a
**percentage of sites that were successfully genotyped** in a given
BovCADD score class, per sample:

This normalization makes burden estimates comparable **within and across
breeds regardless of admixture level**.

---

## 2. BovCADD score classes

| Class            | Threshold      | Interpretation                          |
|------------------|----------------|------------------------------------------|
| Deleterious      | BovCADD > 20   | Top ~1% predicted most damaging variants |
| Benign/Tolerated | BovCADD < 5    | Predicted neutral / low-impact variants  |

For each class, both **homozygosity** and **heterozygosity** percentages
are computed per individual:

| Metric              | Formula                                                        |
|---------------------|------------------------------------------------------------------|
| `Hom_Del_Pct`       | HOM_ALT sites (BovCADD > 20) / genotyped sites (BovCADD > 20) × 100 |
| `Het_Del_Pct`       | HET sites (BovCADD > 20) / genotyped sites (BovCADD > 20) × 100    |
| `Hom_Benign_Pct`    | HOM_ALT sites (BovCADD < 5) / genotyped sites (BovCADD < 5) × 100   |
| `Het_Benign_Pct`    | HET sites (BovCADD < 5) / genotyped sites (BovCADD < 5) × 100      |
| `Burden_Ratio_Hom`  | `Hom_Del_Pct / Hom_Benign_Pct`                                     |
| `Burden_Ratio_Het`  | `Het_Del_Pct / Het_Benign_Pct`                                     |

The **ratio metrics** further isolate whether homozygosity/heterozygosity
is *enriched* at deleterious sites relative to benign background rates —
useful for detecting realized genetic load independent of overall
inbreeding or admixture level.

---

## 3. Repository structure
cattle-genetic-burden/ ├── README.md ├── environment.yml ├── config/ │ ├── config.yaml │ └── sample_breed_map.csv ├── workflow/ │ ├── Snakefile │ └── rules/ │ ├── qc.smk │ ├── bovcadd_annotate.smk │ ├── burden.smk │ └── stats.smk ├── scripts/ │ ├── annotate_bovcadd.sh │ ├── compute_hom_het_burden.py │ ├── merge_chrom_burden.py │ ├── population_stats.py │ └── plot_burden.py ├── resources/ │ └── bovcadd/ # bgzip + tabix-indexed BovCADD score file ├── data/ # input VCFs (gitignored) ├── results/ │ ├── annotated/ │ ├── burden/ │ └── plots/ └── logs/


---

## 4. Installation

```bash
git clone https://github.com/<your-username>/cattle-genetic-burden.git
cd cattle-genetic-burden
conda env create -f environment.yml
conda activate cattle-burden
