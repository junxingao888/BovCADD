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

import io
import sys
import subprocess
from argparse import ArgumentParser

PARSER = ArgumentParser(description=__doc__)
PARSER.add_argument("-i", "--input",
                    help="vcf file to annotate", type=str, required=True)
PARSER.add_argument("-o", "--output",
                    help="Output annotation file, in tsv format",
                    type=str, required=True)
PARSER.add_argument("-a", "--annotation",
                    help="Annotation file, expects a bgzipped, indexed, "
                         "bedgraph (bed.gz)", type=str, required=True)
PARSER.add_argument("-l", "--labels",
                    help="Annotation labels in file, at least 1 should be "
                         "provided, format <1-based column number>=<name>,"
                         "e.g. 4=conservation",
                    type=str, nargs="+")
PARSER.add_argument("-f", "--full",
                    help="prefix for chr in sequence label e.g. chr for "
                         "mm39.chr19 (default None)", action="store_true")
PARSER.add_argument("--overlapping",
                    help="Annotation contains overlapping records, requires "
                         "a different mode in bcftools which is slower and "
                         "memory intensive", action="store_true")
PARSER.add_argument("--logfile", help="Logfile to write coverage information "
                                      "into, (default None)", type=str,
                    default="None")


def build_command(cli_args) -> list:
    """
    Builds BCFtools command used to annotate variants
    :param cli_args: Argparse arguments
    :return: list of command line arguments
    """
    header_lines = []
    annotations = [label.split("=") for label in cli_args.labels]

    columns = ["-"] * (max([int(label[0]) for label in annotations]) - 1)
    columns[0:2] = ["CHROM", "FROM", "TO"]
    for idx, name in annotations:
        idx = int(idx)
        if idx < 3:
            sys.exit(f"The first 3 columns must be Chrom, from/start, end, "
                     f"they cannot be an annotation: {idx} -> {name} !")
        columns[idx - 1] = f"INFO/{name}"
        header_lines.append(f'##INFO=<ID={name},Number=1,Type=String,'
                            f'Description="Score {name}">')

    arguments = ["bcftools", "annotate",
                 "--no-version", "-a", cli_args.annotation,
                 "-c", ",".join(columns), "-O", "v"]
    if not cli_args.overlapping:
        arguments.append("--single-overlaps")
    for header in header_lines:
        arguments.append("-H")
        arguments.append(header)
    arguments.append(cli_args.input)
    return arguments


if __name__ == '__main__':
    args = PARSER.parse_args()
    annotations = [label.split("=") for label in args.labels]
    [sys.exit(f"invalid label: {label}")
     for label in annotations if len(label) != 2]
    labels = [name for idx, name in annotations]

    outfile = open(args.output, "w")

    label_t = "\t".join(labels)
    if args.full:
        outfile.write(f"#Chrom\tPos\tRef\tAlt\t{label_t}\n")
    else:
        outfile.write(label_t + "\n")

    command = build_command(args)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)

    num_var = 0
    num_scores = dict([(label, 0) for label in labels])
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        if line.startswith("#"):
            continue
        num_var += 1
        variant = line.strip().split("\t")

        if args.full:
            outfile.write(f"{variant[0]}\t{variant[1]}\t"
                          f"{variant[3]}\t{variant[4]}\t")

        entries = variant[7].split(";")
        scores = []
        for label in labels:
            found = False
            for entry in entries:
                if label in entry:
                    scores.append(entry.split("=")[1])
                    num_scores[label] += 1
                    found = True
                    break
            if not found:
                scores.append("-")
        outfile.write("\t".join(scores) + "\n")
    outfile.close()

    proc.wait()
    if proc.returncode != 0:
        sys.exit(f"BCFtools returned with a non-zero code: {proc.returncode}")

    if args.logfile != "None":
        logfile = open(args.logfile, "w")
        logfile.write("Ran Command:\n" + " ".join(command) + "\n")
        logfile.write(f"Total variants: {num_var}\nNumber of scores:\n")
        for name, scores in num_scores.items():
            logfile.write(
                f"{name}: {scores} total, coverage: {scores / num_var}\n")
        logfile.close()
