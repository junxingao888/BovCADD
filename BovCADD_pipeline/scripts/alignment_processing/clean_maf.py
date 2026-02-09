#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Date: 04-10-2023

Converts ambiguous nucleotides into N, since maftools does not support the
other IUPAC ambiguous codes. """

import gzip
from argparse import ArgumentParser
# from Bio import AlignIO

IUPAC_AMBIGUOUS = "RYSWKMBDHVryswkmbdhv"

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Input MAF file to filter, can be gzipped",
                    type=str, required=True)
parser.add_argument("-o", "--output",
                    help="Output MAF file to filter, wil be gzipped if "
                         "ending in .gz (default: out.maf)",
                    type=str, default="out.maf")
args = parser.parse_args()

infile = gzip.open(args.input, "rt") \
    if args.input.endswith('.gz') else open(args.input, "r")
outfile = gzip.open(args.output, "wt", compresslevel=1) \
    if args.output.endswith('.gz') else open(args.output, "w")

for line in infile:
    if not line.startswith("s "):
        outfile.write(line)
        continue
    parts = line.split(" ")
    sequence = parts.pop()  # Sequence is last, taken out for filtering
    sequence = sequence.replace('.', 'N')   # <--- ?????,??? .
    for char in IUPAC_AMBIGUOUS:
        sequence = sequence.replace(char, "N")
    outfile.write(" ".join(parts) + " " + sequence)
# alignments = AlignIO.parse(infile, "maf")
# AlignIO.write(alignments, outfile, "maf")
infile.close()
outfile.close()
