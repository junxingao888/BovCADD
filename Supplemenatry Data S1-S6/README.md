# BovCADD: supplementary data

Supplementary data files accompanying:

**Source of genome-wide deleterious variation in a global cattle cohort**

Code and workflows: https://github.com/junxingao888/BovCADD

---

## Overview

BovCADD is a genome-wide deleteriousness score for all 8.1 billion
possible single-nucleotide substitutions in the bovine reference
genome ARS-UCD1.2. Scores are PHRED-scaled by rank, from 1 to 99,
such that values of 10, 20, 30 and 40 correspond to the top 10%, 1%,
0.1% and 0.01% of substitutions genome-wide.

This deposit contains the supplementary data files cited in the
manuscript. Precomputed genome-wide scores are deposited separately
at [DOI / URL].

---

## Files

| File | Contents |
|---|---|
| `Supplemental_Data_S1.xlsx` | ROC-AUC for each annotation individually, by genomic region |
| `Supplemental_Data_S2.xlsx` | Pairwise correlation coefficients among annotations |
| `Supplemental_Data_S3.xlsx` | Feature weights from the L2-regularized logistic regression |
| `Supplemental_Data_S4.xlsx` | BovCADD vs FAETH correlation across score percentile bins |
| `Supplemental_Data_S5.xlsx` | High-burden genes at thresholds of >20, >25, >30, >35 rare deleterious variants per gene |
| `Supplemental_Data_S6.xlsx` | Sample lists: 1000 Bull Genomes, BovCADD training set, OPTIBOV, gene-based burden analyses |

Each file opens with a sheet describing its columns.

---

## Conventions

- **Reference genome:** ARS-UCD1.2 (GCF_002263795.1). Coordinates are 1-based.
- **Scores:** PHRED-scaled BovCADD, range 1-99. Higher values indicate greater predicted deleteriousness.
- **Thresholds used in the manuscript:** high-burden, BovCADD > 20; low-burden, BovCADD < 5.
- **Missing values:** coded as `NA`.
- **Annotation names** follow the schema in Additional file 1: Table S2.

---

## Data sources

Training and analysis data derive from the 1000 Bull Genomes Project
Run 9 (Hayes and Daetwyler 2019) and the LEAP-Agri OPTIBOV project
(Ginja et al. 2025). Access to the underlying sequence data is
governed by those consortia; sample identifiers are listed in
Supplemental Data S6. Annotation sources are listed in Additional
file 1: Table S1.

---

## Licence

[e.g. CC BY 4.0]. If you use these data, please cite the manuscript above.

---

## Contact

Junxin Gao, Wageningen University & Research, [junxin.gao@wur.nl]
