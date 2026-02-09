#!/usr/bin/env python
# -*- coding: ASCII -*-
"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 06-08.2018

This files goes through the ancestral sequence
and compares it with the genome
and the known population variants to identify derived variants.

I make use only of high quality ancestor and assembly
Cases		:1		2		3		4		5
Ancestor	:C		C		C		C		C
Reference	:C		A>0.9	A		A>0.9	A
Alternative	:A>0.9	G		G>0.9	C		A


:Edited by: Seyan Hu
:Date: 19-10-2022
:Edited by: Job van Schipstal
:Date: 27-10-2023
:Usage: See derive_variants.py --help
Changed to better follow PEP8 style and switched to argparse for input.
"""

# Import dependencies
import gzip
import sys
from argparse import ArgumentParser

import pysam

parser = ArgumentParser(description=__doc__)
parser.add_argument('-o', '--output',
                    help='Vcf output file prefix and path, both low '
                         'and high quality variants will be generated '
                         '(default: derived_variants_)',
                    default="derived_variants_")
parser.add_argument('-c', '--chrom',
                    help='investigated Chromosome',
                    required=True)
parser.add_argument('-a', '--ancestor',
                    help='The ancestral sequence file for the chromosome',
                    required=True)
parser.add_argument('-r', '--reference',
                    help='The reference genome file for the chromosome',
                    required=True)
parser.add_argument('-v', '--variants',
                    help='The population variants file for the chormosome',
                    required=True)
parser.add_argument('--only-hq',
                    help='Consider lowercase reference also as '
                         'low_quality output (def include in high_output, '
                         'only look at ancestral quality)',
                    action="store_true")
parser.add_argument('-s', '--start', help='start position of the region to '
                                          'generate variants for',
                    default=0, type=int)

NUCLEOTIDES = ["A", "C", "G", "T"]


def main(args):
    # Define input files for the Ancestral sequence and the reference.
    ancestor_fasta = pysam.Fastafile(args.ancestor)
    ref_fasta = pysam.Fastafile(args.reference)

    # Open population frequency files
    # Note: frequency files are not vcf files
    input_snps = gzip.open(args.variants, "rt") \
        if args.variants.endswith('.gz') else open(args.variants, "r")

    # Check if the ancestor seq is of same length as the reference
    if ancestor_fasta.nreferences != 1:
        sys.exit('There are more than 1 record in the fasta file.')
    ancestor_record = ancestor_fasta.references[0]
    if '.%s' % args.chrom in ancestor_record or \
            'chr%s' % args.chrom in ancestor_record:
        anc_length = ancestor_fasta.get_reference_length(
            ancestor_fasta.references[0])
    else:
        sys.exit(
            'The requested chromosome cannot be found in the record of '
            'the ancestor fasta\n%s\n%s' % ('.%s:' % args.chrom,
                                            ancestor_record))
    if anc_length != ref_fasta.get_reference_length(ref_fasta.references[0]):
        sys.exit(
            'Ancestor fasta and ref fasta have not the same length, '
            'Chromosome: %s' % ('%s' % args.chrom))
    print('Ancestor sequence control: Chr' + str(
        args.chrom) + ', sequence is good')

    # Create output files for upper and lower case variants.
    output_high = open(f"{args.output}_case_upper.vcf", 'w')
    output_low = open(f"{args.output}_case_lower.vcf", 'w')
    for out in [output_low, output_high]:
        out.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL"
                  "\tFILTER\tINFO\n")

    # Split the ancestor and reference sequences in chunks of 350000.
    # Done so that the script reads this in chunks.
    index_outer = (args.start, args.start + 350000)
    input_snps.readline()  # Skip header
    snp_line = input_snps.readline()
    snp_pos = int(snp_line.strip().split('\t')[1]) - 1

    while index_outer[0] <= anc_length:
        # Fetch the reference and ancestral sequence.
        ref_seq = ref_fasta.fetch(ref_fasta.references[0], index_outer[0],
                                  index_outer[1])
        anc_seq = ancestor_fasta.fetch(ancestor_record, index_outer[0],
                                       index_outer[1])

        # Iterate over the ancestral sequence position.
        for i, anc in enumerate(anc_seq):
            ref = ref_seq[i]

            if not args.only_hq:
                ref = ref.upper()

            # If either is lower case, write to low_quality output
            if ref.islower() or anc.islower():
                outfile = output_low
                ref = ref.upper()
                anc = anc.upper()
            else:
                outfile = output_high

            # Replace ref with a more common variant within the population
            # vcf if above threshold. Thus checking if current position in
            # the sequences is equal to the snp position in the frequency
            # files
            if i + index_outer[0] == snp_pos:
                # If we have an snp, this is the reading reference
                # A variant with greater frequency will be taken instead,
                # reference is then ignored
                ref = snp_line.strip().split('\t')[5].split(':')[0].upper()

                # Get next line, if there are still lines
                snp_line = input_snps.readline()
                if snp_line == "":  # No newline character means EOF
                    snp_pos = -1
                else:
                    # -1 to convert to 0-based coordinates
                    snp_pos = int(snp_line.strip().split('\t')[1]) - 1

            # Only look at defined positions in both ref and alt
            if anc not in NUCLEOTIDES or ref not in NUCLEOTIDES:
                continue

            # Write if differing
            if anc != ref:
                outfile.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" % (
                    args.chrom, i + index_outer[0] + 1, ref, anc))

        index_outer = (index_outer[1], index_outer[1] + 350000)

    output_high.close()
    output_low.close()


if __name__ == '__main__':
    main(parser.parse_args())
