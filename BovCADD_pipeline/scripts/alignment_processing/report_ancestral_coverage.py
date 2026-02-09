#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: jgc.vanschipstal@gmail.com
:Date: 03-12-2023
Compares the reference and ancestral genomes and returns coverage statistics.

It is key that an equal amount of ancestor and reference chromosomes are
provided, in the form of individual fasta files.
The order must be the same. e.g.
python <script>.py -a chr1.fa chr2.fa -r ref_chr1.fa ref_chr2.fa
is the intended way to run the script.
"""
# Import dependencies
import sys
from argparse import ArgumentParser
from typing import Generator, Dict

parser = ArgumentParser(description=__doc__)
parser.add_argument("-a", "--ancestor",
                    help="extracted ancestor sequence, 1 or more FASTA files",
                    type=str, required=True, nargs="+")
parser.add_argument("-r", "--reference",
                    help="reference sequence of the species of interest",
                    type=str, required=True, nargs="+")
parser.add_argument("-o", "--outfile", help="Outfile (default stdout)",
                    default="std", type=str)
parser.add_argument("-f", "--fraction",
                    help="File to write just the total fraction covered to "
                         "(Not written by default). Useful to estimate "
                         "number of needed simulations needed later.",
                    type=str, required=False)
args = parser.parse_args()


def yield_sequence(fasta: str) -> Generator[str, None, None]:
    """
    Open fasta file and yield nucleotides in sequence
    :param fasta: str, filename of fasta file
    :return: Yields str of 1 char, individual nucleotides in the fasta file
    """
    with open(fasta, "r") as fasta_f:
        first_l = fasta_f.readline()
        if not first_l.startswith(">"):
            sys.exit(
                f"Expecting fasta files, with a header file starting "
                f"with >, not {first_l}\n")
        for line in fasta_f:
            if line.startswith(">"):
                sys.exit("Only one sequence expected per fasta file!")
            for nt in line.strip():
                yield nt


def compare_sequences(reference: str, ancestor: str) -> Dict[str, int]:
    """
    Compares two sequences and returns alignment statistics
    :param reference: filename and path of reference sequence
    :param ancestor: file name and path of ancestor sequence
    :return: counts for gaps, aligning nt's and mismatches
    """
    count = {"total": 0,
             "gap": 0,
             "ref_undef": 0,
             "low_align": 0,
             "high_align": 0,
             "total_align": 0,
             "high_mismatch": 0,
             "low_mismatch": 0}

    ref_generator = yield_sequence(reference)
    # Iterate through each character in the sequence
    for ref, an in zip(ref_generator, yield_sequence(ancestor)):
        count["total"] += 1
        if ref not in "ACGTacgt" or an not in "ACGTacgt":
            count["gap"] += 1
            if ref not in "ACGTacgt":
                count["ref_undef"] += 1
            continue
        count["total_align"] += 1
        if ref in "ACGT" and an in "ACGT":
            if ref != an:
                count["high_mismatch"] += 1
            count["high_align"] += 1
            continue
        if ref.upper() != an.upper():
            count["low_mismatch"] += 1
        count["low_align"] += 1

    # Verify that end of reference was reached
    if next(ref_generator, None) is not None:
        sys.exit(f"Ancestral sequence does not cover whole reference!\n"
                 f"Generator quit at position {count['total']},\n"
                 f"Ref: {reference}, anc: {ancestor}")
    return count


def print_table(output_dict, file_h, exact=True) -> None:
    """
    Print values in dict of dict as a named table
    :param output_dict: dict of dict, values to print
    :param file_h: openfile-like object to print values to
    :param exact, bool (def T), if false print only 3 digit precision
    :return: None, output is printed to file_h
    """
    first = True
    for (entry, entry_counts) in output_dict.items():
        if first:
            file_h.write("#Name\t" +
                         "\t".join(list(entry_counts.keys())))
            first = False
        file_h.write(f"\n{entry}")
        for (key, value) in entry_counts.items():
            file_h.write(f"\t{value}" if exact else f"\t{value:.3f}")
    file_h.write("\n")


if __name__ == '__main__':
    outfile = sys.stdout if args.outfile == "std" else open(args.outfile, "w")

    if len(args.ancestor) != len(args.reference):
        sys.exit("Error: Equal amount of reference and alt sequences needed!")

    counts = {}
    for ref_f, anc_f in zip(args.reference, args.ancestor):
        anc_n = anc_f.split("/")[-1].split(".")[0]
        sys.stderr.write(f"Calculating coverage for {anc_n}...\n")
        counts[anc_n] = compare_sequences(ref_f, anc_f)

    sys.stderr.write(f"Processing counts...\n")
    totals = {}
    # Sum values in the nested dictionaries to get total counts
    for sub in counts.values():
        for key, sub_count in sub.items():
            totals[key] = sub_count + totals.get(key, 0)
    counts["Total"] = totals

    outfile.write("## Counts by chrom and total:\n")
    print_table(counts, outfile)

    # get fractional values
    counts_frq = {}
    for (chrom, chrom_counts) in counts.items():
        total = float(chrom_counts["total"])
        count_frqs = {}
        for key, sub_count in chrom_counts.items():
            if key != "total":
                count_frqs[key + "_frq"] = sub_count / total
            counts_frq[chrom] = count_frqs

    outfile.write("\n## Frequencies by chrom and total:\n")
    print_table(counts_frq, outfile, False)

    # Write out just fraction covered for use in the pipeline.
    if args.fraction:
        with open(args.fraction, "w") as fraction_f:
            fraction_f.write(f'{counts_frq["Total"]["total_align_frq"]:.4f}')
