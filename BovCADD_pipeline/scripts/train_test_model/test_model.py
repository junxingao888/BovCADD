#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 20-10-2023
:Usage: see test_model.py --help

Loads model chunks from npz files, and the y variable from the metadata csv.
Merges chunks and then scales by standard deviation but not mean centered,
as this would break sparsity. Generates basic statistics for the performance
on the test set of the model that is specified.

Will print the statistics to stdout, progress to stderr.
A roc-auc plot is saved to file if --roc-plot <file> is given.
"""

# Import dependencies
import sys
from argparse import ArgumentParser
from pathlib import Path

from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix, \
    classification_report, roc_curve

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="test set variants, npz file(s), at least one should "
                         "be provided. It is expected that there are "
                         ".meta.csv.gz and columns.csv files present as well",
                    type=str, required=True, nargs="+")
parser.add_argument("--model",
                    help="Pickle file containing the scikit-learn model",
                    type=str, required=True)
parser.add_argument("--scaler",
                    help="Pickle file containing the scikit-learn scaler",
                    type=str, required=True)
parser.add_argument("--roc-plot",
                    help="File to write ROC plot of the model on the test "
                         "set to, (default None, not written). It is "
                         "expected that a format extension (e.g. png or svg) "
                         "supported by pylot is given",
                    type=str, default="None")
parser.add_argument("--columns",
                    help="file (default None, use all columns), "
                         "only train a model based on the "
                         "columns in this file, csv format",
                    type=str, default=None)
parser.add_argument("-n", "--n-jobs",
                    help="int, number of jobs for training the models "
                         "(default 0, as many as needed, "
                         "limited to half the total cores of the device)",
                    type=int, default=0)
parser.add_argument( "--module",
                    help="Location of helper script (optional), "
                         "not needed if in the same folder",
                    type=str, required=False)
args = parser.parse_args()

# Load helper functions.
# Snakemake caching makes this more complicated,
# have to enable importing from another folder
if args.module:
    sys.path.append(Path(args.module).parent)
from data_helper import load_dataset, get_pickle


sys.stderr.write("## Loading model and scaler\n")
model = get_pickle(args.model, LogisticRegression)
scaler = get_pickle(args.scaler, StandardScaler)

# Load dataset
n_jobs = args.n_jobs if args.n_jobs != 0 else len(args.input)
x_test, y_test, columns = load_dataset(args.input,
                                       args.columns,
                                       jobs=args.n_jobs)

sys.stderr.write("## Scaling dataset\n")
x_test = scaler.transform(x_test)
del scaler

# Generate basic statistics on test set
sys.stderr.write("## Calculating scores\n")
print("Confusion matrix:")
y_pred = model.predict(x_test)
print(confusion_matrix(y_test, y_pred))
print("Standard scikit-learn classification report output:")
print(classification_report(y_test, y_pred))

# Calculate ROC-AUC and plot ROC if roc_plot argument was given
y_pred_prob = model.predict_proba(x_test)[::, 1]
fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
auc = roc_auc_score(y_test, y_pred_prob)
print(f"ROC-AUC score: {auc}")
if args.roc_plot != "None":
    plt.plot(fpr, tpr, label="test-set, auc=" + str(auc))
    plt.xlabel("False positive rate (FPR)")
    plt.ylabel("True positive rate (TPR)")
    plt.legend(loc=4)
    plt.savefig(args.roc_plot)
# Save raw ROC data alongside PNG  ### NEW
    out_tsv = Path(args.roc_plot).with_suffix(".tsv")
    import pandas as pd
    pd.DataFrame({
        "fpr": fpr,
        "tpr": tpr
    }).to_csv(out_tsv, sep="\t", index=False)