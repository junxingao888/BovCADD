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
from argparse import ArgumentParser
from typing import Union

import pysam
import numpy as np

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Path to vcf source file (vcf file of the generated "
                         "variants, may be gzipped)", type=str)
parser.add_argument("-r", "--reference", type=str,
                    help="Reference chr file and path with wildcard [CHROM]")
parser.add_argument("-o", "--output",
                    help="Output file (default basic_annotation.tsv)",
                    type=str, default="basic_annotation.tsv")
parser.add_argument("-s", "--dna-shape",
                    help="DNA shape file, def (None -> no shape scores)",
                    type=str, default="None")
parser.add_argument("--include-masked",
                    help="Include isMasked annotation (def False)",
                    action="store_true")
args = parser.parse_args()

# List for transversions and transitions.
TRANSVERSIONS = {('A', 'C'), ('C', 'A'), ('T', 'A'), ('A', 'T'), ('C', 'G'),
                 ('G', 'C'), ('G', 'T'), ('T', 'G')}
TRANSITIONS = {('C', 'T'), ('T', 'C'), ('G', 'A'), ('A', 'G')}

# Size of window to calculate GC and CpG content in.
WINDOWSIZE = 75

HEADER = ['#Chrom', 'Pos', 'Ref', 'Alt', 'isTv', 'GC', 'CpG']
if args.include_masked:
    HEADER.append("isMasked")
if args.dna_shape != "None":
    HEADER.extend(["dnaRoll", "dnaProT", "dnaMGW", "dnaHelT"])

REFERENCE_FILES = {}


## TODO incorporate
# exonTabix = pysam.Tabixfile(options.exons, 'r')


def get_seq_window(chromosome: str, start: int, end: int) -> Union[str, None]:
    """
    Get a sequence window from a specific reference chromosome
    :param chromosome: str, name of chromsome
    :param start: int, start position (0-based)
    :param end: int, end position (0-based, up to not including end)
    :return: str, sequence window, or None if not able to fetch
    """
    ref_tabix, ref_name = REFERENCE_FILES.get(chromosome, (None, None))

    # Load chromosome: Open reference of chr, extract name from index
    if not ref_tabix:
        fasta_f = args.reference.replace("[CHROM]", chromosome)
        ref_tabix = pysam.Fastafile(fasta_f)

        # Read chr number or ncbi identifier from reference index file.
        with open(fasta_f + '.fai', 'r') as ref_idx:
            for line in ref_idx:
                ref_name = line.split('\t')[0]
                break
        REFERENCE_FILES[chromosome] = (ref_tabix, ref_name)
    try:
        return ref_tabix.fetch(ref_name, start, end)
    except:
        pass
    return None


def count_gc_cpg(seq):
    """
    Function for counting GC and CpG sites in a window.
    :param seq: str, window sequence to determine GC and CpG content in.
    :return: tuple of float gc, cpg or -, - if empty sequence
    """
    seq = seq.upper()
    cpg, gc = 0, 0
    count = 0
    lbase = ''
    for pos in range(len(seq)):
        base = seq[pos]
        count += 1
        if base in 'GC':
            gc += 1
        if lbase == 'C' and base == 'G':
            cpg += 1
        lbase = base
    if count > 0:
        return gc / float(count), cpg / (count * 0.5)
    else:
        return '-', '-'


def parse_dna_shape(dnashape_path):
    """
    Build dna shape dictionary based on precalculated scores for all pentamers
    :param dnashape_path: str, file to load from
    :return: dict of dict, dna shape scores by allele and then pentamer
    """
    dna_file = open(dnashape_path, 'r')
    allele_dict_dict = {}
    seq_dict_values = {}
    investigated_allele = ''

    # investigated allele, pentamer sequence, Roll, Prot, MGW, HelT
    for lines in dna_file:
        lines = lines.strip()
        if not lines:
            continue
        if investigated_allele == '':
            investigated_allele = lines[0]

        if investigated_allele == lines[0]:
            line = lines.strip().split('\t')
            seq_dict_values[line[1]] = line[2:]
        else:
            allele_dict_dict[investigated_allele] = seq_dict_values
            seq_dict_values = {}
            line = lines.strip().split('\t')
            investigated_allele = line[0]
            seq_dict_values[line[1]] = line[2:]
    allele_dict_dict[investigated_allele] = seq_dict_values
    return allele_dict_dict


def get_shape_score(chrom, start, ref_allele, alt_allele, parsed_shape):
    diff = abs(len(ref_allele) - len(alt_allele))

    seq_ref = get_seq_window(chrom, start - 3, start + 2 + diff)
    if seq_ref is None:
        print(f"[WARNING] get_seq_window returned None for {chrom}:{start-3}-{start+2+diff}")
        return None
    seq_ref = seq_ref.upper()

    if len(ref_allele) > len(alt_allele):
        new_allele_ref = ref_allele
        new_allele_alt = alt_allele + seq_ref[3:-2]
        begin_seq_ref = seq_ref[:2]
        end_seq_ref = seq_ref[-2 - diff:-diff]
        begin_seq_alt = seq_ref[:2]
        end_seq_alt = seq_ref[-2:]

    elif len(ref_allele) < len(alt_allele):
        new_allele_ref = ref_allele + seq_ref[3:-2]
        new_allele_alt = alt_allele
        begin_seq_ref = seq_ref[:2]
        end_seq_ref = seq_ref[-2:]
        begin_seq_alt = seq_ref[:2]
        end_seq_alt = seq_ref[-2 - diff:-diff]

    else:
        new_allele_ref = ref_allele
        new_allele_alt = alt_allele

    # Score indel
    if (len(ref_allele) != 1) or (len(alt_allele) != 1) or ('-' in ref_allele) or ('-' in alt_allele):
        if ('N' in new_allele_alt) or ('N' in new_allele_ref):
            return None
        query_ref = []
        query_alt = []
        shape_scores = []

        nucleotides = ["", "", "", "", ""]
        for i, nuc in enumerate(new_allele_ref):
            switch = 0
            nucleotides[2] = nuc
            try:
                nucleotides[3] = new_allele_ref[i + 1]
                switch = 1
            except IndexError:
                nucleotides[3] = end_seq_ref[0]
            try:
                nucleotides[4] = new_allele_ref[i + 2]
            except IndexError:
                nucleotides[4] = end_seq_ref[1 - switch]
            query = begin_seq_ref[i:] + "".join(nucleotides)
            nucleotides[0] = nucleotides[1]
            nucleotides[1] = nucleotides[2]
            query_ref.append(query)

        shape_ref = []
        for string in query_ref:
            try:
                shape_ref.append(parsed_shape[string[2]][string])
            except KeyError:
                return None

        nucleotides = ["", "", "", "", ""]
        for i, nuc in enumerate(new_allele_alt):
            switch = 0
            nucleotides[2] = nuc
            try:
                nucleotides[3] = new_allele_alt[i + 1]
                switch = 1
            except IndexError:
                nucleotides[3] = end_seq_alt[0]
            try:
                nucleotides[4] = new_allele_alt[i + 2]
            except IndexError:
                nucleotides[4] = end_seq_alt[1 - switch]
            query = begin_seq_alt[i:] + "".join(nucleotides)
            nucleotides[0] = nucleotides[1]
            nucleotides[1] = nucleotides[2]
            query_alt.append(query)

        shape_alt = []
        for string in query_alt:
            try:
                shape_alt.append(parsed_shape[string[2]][string])
            except KeyError:
                return None

        if len(shape_alt) == len(shape_ref):
            for aa, bb in zip(shape_alt, shape_ref):
                shape_scores.append([float(a) - float(b) for a, b in zip(aa, bb)])
            return np.apply_along_axis(np.mean, 0, np.array(shape_scores))
        else:
            return None

    # Score SNV
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        seq_ref = get_seq_window(chrom, start - 3, start + 2).upper()
        if 'N' in seq_ref or len(seq_ref) != 5:
            return None
        if seq_ref[2] != ref_allele:
            return None
        try:
            shape_ref = parsed_shape[ref_allele][seq_ref]
            seq_alt = seq_ref[:2] + alt_allele + seq_ref[3:]
            shape_alt = parsed_shape[alt_allele][seq_alt]
            return [float(a) - float(b) for a, b in zip(shape_alt, shape_ref)]
        except KeyError:
            return None

    return None


if __name__ == '__main__':
    outfile = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "w")
    outfile.write("\t".join(HEADER))

    infile = gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "r")

    dna_shape = parse_dna_shape(args.dna_shape) if args.dna_shape != "None" else None

    # Annotate each variant
    for line in infile:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        elems = line.split("\t")
        chrom = elems[0]
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        pos = int(elems[1])
        ref = elems[3]
        alt = elems[4]

        is_tv = (ref, alt) in TRANSVERSIONS
        outfile.write(f"\n{chrom}\t{pos}\t{ref}\t{alt}\t{is_tv}")

        sequence = get_seq_window(chrom, pos - WINDOWSIZE - 1, pos + WINDOWSIZE)
        if sequence:
            is_masked = 0 if sequence[75].isupper() else 1
            gc, cpg = count_gc_cpg(sequence)
            outfile.write(f"\t{gc:.3f}\t{cpg:.4f}")
        else:
            is_masked = "-"
            outfile.write("\t-\t-")

        if args.include_masked:
            outfile.write(f"\t{is_masked}")

        if dna_shape:
            scores = get_shape_score(chrom, pos, ref, alt, dna_shape)
            if scores:
                for score in scores:
                    outfile.write(f"\t{score:.3f}")
            else:
                outfile.write("\t-\t-\t-\t-")

    outfile.write("\n")
    outfile.close()
    infile.close()
