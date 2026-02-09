#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 15-12-2023
:Usage: see basic_annotation.py --help

Generate basic annotation columns.
These can be derived from just the vcf and the reference genome.

These functions were copied (with minor modifications) from Gross et al.:
From pCADD VEP processing script:
- count_gc_cpg
From mCADD "Mouse-feature_annotate.py"
- parse_dna_shape
- get_shape_score (previously part of a larger function)
"""

# Import dependencies
import gzip
import math
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas
import pandas as pd
import pysam
import numpy as np

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Path to input .tsv file containing all annotations",
                    type=str, nargs="+")
parser.add_argument("-o", "--output",
                    help="Output image file",
                    type=str, default="cadd_bins.png")
args = parser.parse_args()

if __name__ == '__main__':
    ABBREVIATIONS = ["SG", "CS", "NS", "SN", "SL", "S", "U5", "U3", "R",
                     "IG", "NC", "I", "UP", "DN", "O"]
    COLUMNS = ["Bin", "Stop gained", "Canonical splice", "Non-synonymous",
               "Synonymous", "Stop lost", "Splice site", "5' UTR", "3' UTR",
               "Regulatory", "Intergenic", "Noncoding change", "Intronic",
               "Upstream", "Downstream", "Unknown"]

    bins = []
    for i in range(52):
        bins.append([0] * len(ABBREVIATIONS))

    # Combine all counts into a single matrix
    for infile in args.input:
        with open(infile, "r") as count_file:
            for line in count_file:
                bin_id, scores = line.strip().split("\t")
                bin_id = int(bin_id)
                scores = scores.split(",")
                for idx, score in enumerate(scores):
                    bins[bin_id][idx] += int(score)

    # Normalise counts


    # convert count to proportion
    df_bins = []
    for bin, counts in enumerate(bins):
        total = sum(counts)
        if total == 0:
            total = 1 # TODO remove
        fractions = [bin]
        for idx in range(len(counts)):
            fractions.append((counts[idx] + 0.07) / total) # TODO remove + 1
        df_bins.append(fractions)

    df = pd.DataFrame(df_bins, columns=COLUMNS)
    print(df)
    df.plot(x='Bin', kind='bar', stacked=True, figsize=(9, 3))

    plt.savefig(args.output)