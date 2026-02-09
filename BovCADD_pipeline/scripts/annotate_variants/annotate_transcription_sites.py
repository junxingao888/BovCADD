#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 05-02-2024

The code for determining the minimum transcription site start and end distances
was taken from Christian Gross' Mouse-feature_annotate.py and moved
into this stand-alone script.
It finds the minimum distance to a transcription start and end site, with
a maximum of 10Mb.

:Usage: see annotate_transcription_sites.py --help
"""

import sys
import os
import gzip
from argparse import ArgumentParser

import pysam

PARSER = ArgumentParser(description=__doc__)
PARSER.add_argument("-i", "--input",
                    help="vcf file to annotate", type=str, required=True)
PARSER.add_argument("-o", "--output",
                    help="Output annotation file, in tsv format",
                    type=str, required=True)
PARSER.add_argument("-a", "--annotation",
                    help="Annotation file, expects a bgzipped, indexed, "
                         "preprocessed gtf file (gtf.gz)",
                    type=str, required=True)
PARSER.add_argument("-f", "--full",
                    help="prefix for chr in sequence label e.g. chr for "
                         "mm39.chr19 (default None)", action="store_true")
PARSER.add_argument("-p", "--prefix",
                    help="prefix for chr in sequence label e.g. chr for "
                         "chr19 (default no prefix)", default="")

if __name__ == '__main__':
    args = PARSER.parse_args()
    outfile = open(args.output, "w")

    # open gtf file, allowing for seeking via the index
    if (not os.path.exists(args.annotation) or
            not os.path.exists(args.annotation + ".tbi")):
        sys.exit("vcf source file: Require valid path to compressed tabix "
                 "file and tabix index file.")
    gtf_tabix = pysam.Tabixfile(args.annotation, 'r')

    # write header
    if args.full:
        outfile.write("#Chrom\tPos\tRef\tAlt\t")
    outfile.write("minDistTSS\tminDistTSE\n")

    maxvalDIST = 100000

    # loop through all variants
    variant_file = gzip.open(args.input, "rt")
    for line in variant_file:
        if line.startswith("#"):
            continue
        variant = line.strip().split("\t")
        chrom = variant[0]
        pos = int(variant[1])

        if args.full:
            outfile.write(f"{chrom}\t{pos}\t"
                          f"{variant[3]}\t{variant[4]}\t")

        minDistTSS, minDistTSE = None, None
        for cmaxval in (maxvalDIST / 100,
                        maxvalDIST / 10,
                        maxvalDIST):
            for elem in gtf_tabix.fetch(args.prefix + chrom,
                                        max(0, pos - cmaxval - 1),
                                        pos + cmaxval + 1):
                elem = elem.rstrip("\n").split('\t')
                start, end = int(elem[3]), int(elem[4])
                if "-" in elem[3]:
                    start, end = end, start
                val = min(max(start - pos + 1, pos - start), maxvalDIST)
                if minDistTSS is None or val < minDistTSS:
                    minDistTSS = val
                val = min(max(end - pos + 1, pos - end), maxvalDIST)
                if minDistTSE is None or val < minDistTSE:
                    minDistTSE = val
            if minDistTSS and minDistTSE:
                outfile.write(f"{minDistTSS}\t{minDistTSE}\n")
                break
        if not minDistTSS or not minDistTSE:
            outfile.write(f"{maxvalDIST}\t{maxvalDIST}\n")

    outfile.close()
    variant_file.close()
