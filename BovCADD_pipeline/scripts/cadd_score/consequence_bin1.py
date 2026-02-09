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

import gzip
import math
import sys
from argparse import ArgumentParser

import pandas as pd
import pysam

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Path to input .tsv file containing all annotations",
                    type=str, required=True)
parser.add_argument("-o", "--output",
                    help="Output file (default consequence_bins.tsv)",
                    type=str, default="basic_annotation.tsv")
parser.add_argument("-a", "--annotation",
                    help="CADD score annotation file, bed.gz format",
                    type=str, required=True)
args = parser.parse_args()

if __name__ == '__main__':
    # Output file handler (gzip or text)
    outfile = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "w")

    # Load annotation file and ensure correct types/format
    df = pd.read_csv(args.input, sep='\t', na_values=['-'])

    df["#Chrom"] = df["#Chrom"].astype(str).str.strip()
    df["Ref"] = df["Ref"].astype(str).str.strip().str.upper()
    df["Alt"] = df["Alt"].astype(str).str.strip().str.upper()
    df["Pos"] = df["Pos"].astype(int)

    # Open CADD file with tabix
    cadd_tabix = pysam.TabixFile(args.annotation, 'r')

    ABBREVIATIONS = ["SG", "CS", "NS", "SN", "SL", "S", "U5", "U3", "R",
                     "IG", "NC", "I", "UP", "DN", "O"]

    bins = []
    for i in range(52):
        bins.append([0] * len(ABBREVIATIONS))

    for index, row in df.iterrows():
        chrom = str(row["#Chrom"]).strip()
        pos = int(row["Pos"])
        ref = str(row["Ref"]).strip().upper()
        alt = str(row["Alt"]).strip().upper()
        found = False
        for elem in cadd_tabix.fetch(chrom, pos - 1, pos):
            elem = elem.rstrip("\n").split('\t')
            # Make sure CADD file's REF/ALT are also formatted properly
            cadd_pos = int(elem[1])
            cadd_ref = elem[2].strip().upper()
            cadd_alt = elem[3].strip().upper()
            if cadd_pos == pos and ref == cadd_ref and alt == cadd_alt:
                score = float(elem[5])
                consequence = row["Consequence"]
                bin_index = min(math.floor(score), 51)  # cap at 51
                bins[bin_index][ABBREVIATIONS.index(consequence)] += 1
                found = True
                break
        if not found:
            print(f"No score found for: {chrom}, {pos}, {ref}, {alt}", file=sys.stderr)
            # Optionally: sys.exit(1) if you want to halt, otherwise continue

    # Write output
    for bin_index, counts in enumerate(bins):
        outfile.write(str(bin_index) + "\t" +
                      ",".join([str(count) for count in counts]) + "\n")

    outfile.close()
