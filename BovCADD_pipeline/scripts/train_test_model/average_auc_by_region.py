#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 5-12-2023
:Usage: see <script>.py --help

Reproduces the prediction performance plot of (only) CADD as in the papers
by Christian Gross. The test dataset is loaded from file, and the predicted
probabilities are loaded from a pre-trained model. Subsets of
the data are taken representing various genomic subsets and for each a
ROC-AUC score is calculated, these are finally plotted in a pyplot.
"""

# Import dependencies
import sys
from pathlib import Path
from argparse import ArgumentParser

import pandas
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import roc_auc_score

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Input npz file(s), testset to generate plot from, "
                         "one for each model. At least one needed",
                    type=str, required=True, nargs="+")
parser.add_argument("-p", "--predictions",
                    help="CSV file(s), containing each model's "
                         "predicted probabilities. must match the input.",
                    type=str, required=True, nargs="+")
parser.add_argument("-o", "--output-plot",
                    help="File to write plot of the model on the test "
                         "set to. It is expected that a format extension "
                         "(e.g. png or svg) supported by pylot is given",
                    type=str, required=True)
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


def get_score(dataframe):
    """
    Calculates the roc_auc score from the values in a dataframe
    :param dataframe: Pandas datframe to extract values from
    :return: float, ROC-AUC score
    """
    return roc_auc_score(dataframe["y_true"], dataframe["y_pred_prob"])


if len(args.input) != len(args.predictions):
    sys.exit(f"Equal amount of datafiles and predictions required!")

desired = ["IND_cDNApos", "IND_CDSpos", "IND_SIFTval"]
auc_scores = {
    "I": [], "II": [], "III": [], "IV": [], "V": [], "VI": [], "VII": [],
}
for infile, predictions in zip(args.input, args.predictions):
    # Get needed columns, these are indicator columns for missing values
    data_matrix, meta_data, cols = load_npz_with_meta([infile], desired)
    del meta_data
    df = pandas.DataFrame(data_matrix.todense(), columns=cols)
    del data_matrix, cols

    # Load predictions
    preds = pandas.read_csv(predictions, sep=',', header=0)
    df = pandas.concat([df, preds], axis=1)
    del preds

    sys.stderr.write(f"## Calculating ROC-AUC scores\n")
    df_cdna = df.loc[df['IND_cDNApos'] == 0.0]
    df_cds = df.loc[df['IND_CDSpos'] == 0.0]
    for key, value in [
        ["I", get_score(df)],
        ["II", get_score(df.loc[df['IND_cDNApos'] != 0.0])],
        ["III", get_score(df_cdna)],
        ["IV", get_score(df_cdna.loc[df_cdna['IND_CDSpos'] != 0.0])],
        ["V", get_score(df_cds)],
        ["VI", get_score(df_cds.loc[df_cds['IND_SIFTval'] != 0.0])],
        ["VII", get_score(df.loc[df['IND_SIFTval'] == 0.0])]
    ]:
        auc_scores[key].append(value)

# Get mean and standard deviation
auc_means = {}
auc_error = {}
for key in auc_scores.keys():
    auc_means[key] = np.mean(auc_scores[key])
    auc_error[key] = np.std(auc_scores[key])

names = ["Total", "non-CDNA", "cDNA", "non-cds",
         "cds", "synonymous", "missense"]
for i, value in enumerate(auc_means.items()):
    print(f"{value[0]} {names[i]}: {value[1]:.3f}±{auc_error[value[0]]:.3f}: "
          f"{auc_scores[value[0]]}")

sys.stderr.write(f"## plotting scores\n")
x = list(auc_means.keys())
y = list(auc_means.values())
err = list(auc_error.values())

plt.figure(figsize=(6, 4))

colors = ["#003780", "#072CA8", "#171DFF", "#3A3FFF", "#5D62FE", "#8084FE",
          "#A3A6FE"]
# colors = ["#0E9ABE", "#10ABD2", "#13BBE7", "#15CCFB", "#34E0FF", "#58F1FF",
#          "#85FCFF"]
plt.bar(x, y, color=colors, width=6 / len(x))
plt.errorbar(x, y, yerr=err, fmt="none", color="#000000", elinewidth=4)

# add value labels
for i in range(len(x)):
    plt.text(i, y[i] + err[i] + 0.01, f"{y[i]:.3f}\n±{err[i]:.3f}",
             ha='center')

plt.title("ROC-AUC scores by genomic region")

plt.xlabel("Genomic region")
plt.ylabel("ROC-AUC")
plt.ylim(0.5, 1.1)

plt.savefig(args.output_plot)
