#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 05-10-2023

Wrapper for annotating variants using BCFTools
Post-processes the annotated VCF into a .tsv file with the scores.

CAUTION: The script assumes a sorted and indexed bedgraph annotation file is
provided. BCFtools is ran with --single-overlaps, not taking overlapping
records in the annotation file into full consideration. This can be enabled
by running the script with --full but this mode requires notably more time
and memory.

A 30GB (uncompressed size) BED file can be used to annotate a million
variants within 30 seconds, using MB's of memory. Doing so with without the
--single-overlaps in BCFtools took 5 minutes and used up to 65GB of memory.

:Usage: see bed_annotation_wrapper.py --help
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
PARSER.add_argument("-b", "--bed",
                    help="Annotation file, expects a bgzipped, indexed, "
                         "bedgraph (bed.gz)", type=str, required=True)
PARSER.add_argument("-a", "--annotations",
                    help="Annotations to extract from the file, "
                         "at least 1 should be provided, separated by spaces "
                         "format: <name>=<type>:<bed-column>. "
                         "Name will be used as the header for the annotation "
                         "type is unique, count, min or max operator",
                    type=str, nargs="+")
PARSER.add_argument("-f", "--full",
                    help="prefix for chr in sequence label e.g. chr for "
                         "mm39.chr19 (default None)", action="store_true")
PARSER.add_argument("--zero-based",
                    help="The annotation file is zero-based half open "
                         "instead of 1-based coordinates.",
                    action="store_true")

if __name__ == '__main__':
    args = PARSER.parse_args()
    annotations = [annotation.split("=") for annotation in args.annotations]
    names = [annotation[0] for annotation in annotations]
    types = [annotation[1].split(":") for annotation in annotations]
    outfile = open(args.output, "w")

    offset = 1 if args.zero_based else 0

    # open bed file, allowing for seeking via the index
    if not os.path.exists(args.bed) or not os.path.exists(args.bed + ".tbi"):
        sys.exit("vcf source file: Require valid path to compressed tabix "
                 "file and tabix index file.")
    bed_tabix = pysam.Tabixfile(args.bed, 'r')

    # write header
    label_t = "\t".join(names)
    if args.full:
        outfile.write("#Chrom\tPos\tRef\tAlt\t")
    outfile.write(label_t + "\n")

    # loop through all variants
    variant_file = gzip.open(args.input, "rt")
    for line in variant_file:
        line = line.strip()
        if line.startswith("#"):
            continue
        variant = line.split("\t")
        chrom = variant[0]
        pos = int(variant[1])

        if args.full:
            outfile.write(f"{chrom}\t{pos}\t"
                          f"{variant[3]}\t{variant[4]}\t")

        # search for matches:
        start = pos - 1 - offset
        end = pos - offset

        """print(f"fetching for {chrom}: {start}-{end}")
        for output in bed_tabix.fetch(chrom, start, end):
            print(output)
            out = output.split('\t')
            if (int(out[1]) > start) & (int(out[1]) <= end):
                matches.append(out)
        """
        matches = [match.split("\t") for
                   match in bed_tabix.fetch(chrom, start, end)]
        # print(f"Fetched {len(matches)} for {chrom}: {start}-{end}")
        if len(matches) == 0:
            outline = "-\t" * len(types)
            outfile.write(outline[:-1] + "\n")
            continue

        scores = []
        for annotation_type, index in types:
            index = int(index)
            values = [match[index - 1] for match in matches]
            if annotation_type == "unique":
                scores.append(len(set(values)))
            elif annotation_type == "count":
                scores.append(len(values))
            else:
                values = [float(value) for value in values]
                if annotation_type == "max":
                    scores.append(max(values))
                elif annotation_type == "min":
                    scores.append(min(values))
                else:
                    sys.exit(f"Unknown annotation type: {annotation_type}")

        outline = ""
        for score in scores:
            if isinstance(score, str):
                outline += score + "\t"
            else:
                outline += format(score, ".4g") + "\t"
        outfile.write(outline[:-1] + "\n")
    outfile.close()
    variant_file.close()
