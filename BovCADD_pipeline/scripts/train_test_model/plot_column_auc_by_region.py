#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 18-10-2023
:Usage: see <script>.py --help

Reproduces the prediction performance plot of (only) CADD as in the papers
by Christian Gross. The test dataset is loaded from file, and the predicted
probabilities are obtained using the pretrained scaler and model. Subsets of
the data are taken representing various genomic subsets and for each a
ROC-AUC score is calculated, these are finally plotted in a pyplot.
"""

# Import dependencies
import sys
from pathlib import Path
from argparse import ArgumentParser

import pandas
from matplotlib import pyplot as plt
from sklearn.metrics import roc_auc_score

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Input npz, testset to generate plot from",
                    nargs="+", required=True)
parser.add_argument("-c", "--columns",
                    help="columns to generate roc-auc plots for.",
                    nargs="+", required=True)
parser.add_argument("-o", "--output-plot",
                    help="File to write plot of the model on the test "
                         "set to. It is expected that a format extension "
                         "(e.g. png or svg) supported by pylot is given.,"
                         "Multiple files should be given if more than 1"
                         "column should be tested",
                    nargs="+", required=True)
parser.add_argument("-m", "--module",
                    help="Location of helper script (optional), "
                         "not needed if in the same folder",
                    type=str, required=False)
args = parser.parse_args()

# Load helper functions.
# Snakemake caching makes this more complicated,
# have to enable importing from another folder
if args.module:
    sys.path.append(Path(args.module).parent)
from data_helper import load_npz_with_meta


def get_score(score_col, dataframe):
    score_raw = dataframe[score_col]
    score = (score_raw - score_raw.min()) / (score_raw.max() - score_raw.min())
    
    df_valid = dataframe[["y_true"]].copy()
    df_valid["score"] = score
    df_valid = df_valid.dropna()

    if df_valid.empty:
        return float("nan")

    return roc_auc_score(df_valid["y_true"], df_valid["score"])

def get_counts(dataframe):
    """
    Calculates the roc_auc score from the values in a dataframe
    :param dataframe: Pandas dataframe to extract values from
    :return: tuple of float, count for y0 and y1
    """
    df_sim = dataframe.loc[dataframe['y_true'] != 0.0]
    df_der = dataframe.loc[dataframe['y_true'] == 0.0]
    return len(df_sim), len(df_der)

desired = ["IND_cDNApos", "IND_CDSpos", "IND_SIFTval"]
desired.extend(args.columns)

# Load data, retrieve y from metadata csv
data, meta_data, cols = load_npz_with_meta(args.input, desired=desired)
y_true = meta_data[["y"]].to_numpy(dtype="float").flatten()
del meta_data
df = pandas.DataFrame(data.todense(), columns=cols)
del data
df["y_true"] = y_true

sys.stderr.write(f"## Calculating ROC-AUC scores\n")
df_cdna = df.loc[df['IND_cDNApos'] == 0.0]
df_cds = df.loc[df['IND_CDSpos'] == 0.0]

names = ["Total", "non-CDNA", "cDNA", "non-cds",
         "cds", "synonymous", "missense"]

sim_count, der_count = get_counts(df)
total_count = sim_count + der_count
print("Variant counts by type:")
print("ID\tName\tsimulated_frq\tderived_frq\ttotal_frq"
      "\tsim_count\tder_count\ttotal_count")
counts = {
        "I": get_counts(df),
        "II": get_counts(df.loc[df['IND_cDNApos'] != 0.0]),
        "III": get_counts(df_cdna),
        "IV": get_counts(df_cdna.loc[df_cdna['IND_CDSpos'] != 0.0]),
        "V": get_counts(df_cds),
        "VI": get_counts(df_cds.loc[df_cds['IND_SIFTval'] != 0.0]),
        "VII": get_counts(df.loc[df['IND_SIFTval'] == 0.0])
}
for i, (name, count) in enumerate(counts.items()):
    combined = count[0] + count[1]
    print(f"{name}\t{names[i]}\t"
          f"{count[0]/sim_count:.3f}\t{count[1]/der_count:.3f}\t"
          f"{combined/total_count:.3f}\t"
          f"{count[0]}\t{count[1]}\t{combined}")

for file, col in zip(args.output_plot, args.columns):
    print(f"Generating ROC-AUC scores for column {col}")
    auc_scores = {
        "I": get_score(col, df),
        "II": get_score(col, df.loc[df['IND_cDNApos'] != 0.0]),
        "III": get_score(col, df_cdna),
        "IV": get_score(col, df_cdna.loc[df_cdna['IND_CDSpos'] != 0.0]),
        "V": get_score(col, df_cds),
        "VI": get_score(col, df_cds.loc[df_cds['IND_SIFTval'] != 0.0]),
        "VII": get_score(col, df.loc[df['IND_SIFTval'] == 0.0]),
    }
    for i, value in enumerate(auc_scores.items()):
        print(f"{value[0]}\t{names[i]}\t{value[1]:.3f}")

    sys.stderr.write(f"## plotting scores for {col}\n")
    x = list(auc_scores.keys())
    y = list(auc_scores.values())

    plt.figure(figsize=(6, 4))

    colors = ["#002452", "#072CA8", "#171DFF", "#3A3FFF", "#5D62FE", "#8084FE",
              "#A3A6FE"]
    plt.bar(x, y, color=colors, width=6/len(x))

    # add value labels
    for i in range(len(x)):
        plt.text(i, y[i] + 0.01, f"{y[i]:.3f}", ha='center')

    plt.title(f"ROC-AUC scores by genomic region for {col}")

    plt.xlabel("Genomic region")
    plt.ylabel("ROC-AUC")
    plt.ylim(0.5, 1.0)

    plt.savefig(file)
