# BovCADD

**Source of genome-wide deleterious variation in a global cattle cohort**

Junxin Gao\*, Martijn F.L. Derks, Job G.C. van Schipstal, Ying Liu,
Etske Bijl, Catarina Ginja, Juha Kantanen, Nasser Ghanem,
Donald Kugonza, Mahlako Makgahlela, Martien A.M. Groenen,
Henk Bovenhuis, Richard P.M.A. Crooijmans\*

\*Corresponding authors

Preprint: https://www.biorxiv.org/content/10.64898/2026.09.07.749866v1

---

## Abstract

Identifying deleterious DNA changes underpins efforts to improve animal
health, welfare and sustainable breeding. In cattle, current variant
prioritization focuses on coding changes, uses single annotation types
and gives limited resolution in non-coding sequence.

We developed BovCADD (bovine Combined Annotation-Dependent Depletion),
a nucleotide-level deleteriousness score for substitutions in
*Bos taurus* and *Bos indicus*, combining evolutionary constraint,
sequence context, epigenetic and regulatory annotations, and gene and
protein features. A logistic regression model trained on 41.9 million
high-frequency derived alleles from about 3,700 cattle, contrasted with
context-matched simulated variants, scored all 8.1 billion possible
substitutions. BovCADD distinguished known pathogenic variants from
background variation, discriminated among variants within the same
consequence class, and scored intronic and intergenic sites.
Aggregating scores identified genes carrying rare deleterious variation
and revealed elevated genetic load at trait-relevant loci and in
bottlenecked, intensively selected populations. BovCADD provides the
first genome-wide, nucleotide-resolution measure of deleteriousness in
cattle, extending variant interpretation into non-coding sequence and
linking variant-level prioritization with population-level patterns of
mutational burden.

---

## Precomputed scores

Scores for all 8.1 billion possible single-nucleotide substitutions are
available for the ARS-UCD1.2 and ARS-UCD2.0 reference assemblies:

https://figshare.com/s/732b5e4306797bd66462

Scores are PHRED-scaled by rank, from 1 to 99. Values of 10, 20, 30 and
40 correspond to the top 10%, 1%, 0.1% and 0.01% of substitutions
genome-wide. Thresholds used in the manuscript: high-burden,
BovCADD > 20; low-burden, BovCADD < 5.

---

## Repository contents

| Path | Contents |
|---|---|
| `annotation/` | Annotation extraction and harmonisation into the 88-feature schema |
| `training/` | Derived-allele calling, variant simulation, model training and cross-validation |
| `scoring/` | Genome-wide scoring and PHRED scaling |
| `analysis/` | OMIA validation, gene-based burden, genetic load, FAETH comparison |
| `figures/` | R scripts reproducing manuscript figures |

---

## Supplementary data

| File | Contents |
|---|---|
| Supplemental Data S1 | ROC-AUC for each annotation individually, by genomic region |
| Supplemental Data S2 | Pairwise correlation coefficients among annotations |
| Supplemental Data S3 | Feature weights from the L2-regularized logistic regression |
| Supplemental Data S4 | BovCADD vs FAETH correlation across score percentile bins |
| Supplemental Data S5 | High-burden genes at thresholds of >20, >25, >30, >35 rare deleterious variants per gene |
| Supplemental Data S6 | Sample lists: 1000 Bull Genomes, BovCADD training set, OPTIBOV, gene-based burden analyses |

---

## Data sources

Training and analysis data derive from the 1000 Bull Genomes Project
Run 9 and the LEAP-Agri OPTIBOV project. Access to the underlying
sequence data is governed by those consortia; sample identifiers are
listed in Supplemental Data S6.

FAETH scores mapped to ARS-UCD1.2 were obtained from Zenodo
(https://doi.org/10.5281/zenodo.20667548).

---

## Citation

If you use BovCADD, please cite the preprint above.

---

## Contact

Junxin Gao, junxin.gao@wur.nl
Richard Crooijmans, richard.crooijmans@wur.nl
Wageningen University & Research
