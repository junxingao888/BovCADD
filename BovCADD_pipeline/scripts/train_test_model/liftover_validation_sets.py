# Import dependencies
import sys
from argparse import ArgumentParser

import pandas
import numpy as np

# Create and process cli
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="input variant csv, as obtained from the mCADD paper",
                    type=str, required=True)
parser.add_argument("-l", "--lifted",
                    help="output of the hgLiftOver tool, (default None, "
                         "will generate the input for liftover)",
                    type=str, default="None")
parser.add_argument("-k", "--known-vcf",
                    help="output vcf file location, for known variants",
                    type=str, default="None")
parser.add_argument("-c", "--control-vcf",
                    help="output vcf file location, for the control variants",
                    type=str, default="None")
args = parser.parse_args()

NUCLEOTIDES = ["A", "C", "T", "G"]


def decode(n_tuple, prefix):
    """
    Retrieves nucleotide from one-hot-encoded nt columns.
    :param n_tuple: named tuple, row to find column in
    :param prefix: str, prefix of column name
    :return: str nucleotide or numpy's Nan if none found
    """
    for nt in NUCLEOTIDES:
        if getattr(n_tuple, prefix + nt) != 0:
            return nt
    sys.stderr.write(f"Unable to find nt in {tuple}\n")
    return np.NaN


# Open DF in pandas.
df = pandas.read_csv(args.input, sep=',', na_values=['-'])
print(df.head())
df = df.rename(
    columns={"#Chrom": "Chrom"})  # The # is a problem for named tuples

if args.lifted == "None":
    with open(args.input + ".pos.txt", "w") as out_f:
        for row in df.itertuples(index=False):
            pos = getattr(row, "Pos")
            out_f.write(f"chr{getattr(row, 'Chrom')}:{pos}-{pos}\n")
        sys.stderr.write(f"No liftover output provided using the --lifted "
                         f"parameters, therefore the input positions boxes "
                         f"for hgLiftOver were written to the file "
                         f"{args.input}.pos.txt\n")
        sys.stderr.write(f"Number of variants: {len(df.index)}\n")
        sys.exit()

# We have lifted so load them
new_pos = []
with open(args.lifted) as in_f:
    for line in in_f:
        chrom_start = line.split("-")[0].split(":")
        if not len(chrom_start) == 2:
            sys.exit(f"Found line without chr and pos {line}")
        new_pos.append(chrom_start)

if len(new_pos) != len(df.index):
    sys.exit(f"Difference in lenght between dataset and liftover,"
             f"did not all variants get lifted over?")
# TODO could allow filtering out of un-lifted variants, not needed currently

control_variants = []
known_variants = []
# Go through file row by row
idx = -1
for row in df.itertuples(index=False):
    idx += 1
    # get new position
    chrom, pos = new_pos[idx]

    # Decode one-hot-encoded ref and alt, store in dict
    variant = {"chrom": chrom,
               "pos": pos,
               "ref": decode(row, "Ref_"),
               "alt": decode(row, "Alt_")}
    if getattr(row, "y") == 0:
        control_variants.append(variant)
    else:
        known_variants.append(variant)

for file, variants in zip((args.known_vcf, args.control_vcf),
                          (known_variants, control_variants)):
    if file != "None":
        with open(file, "w") as out_f:
            out_f.write("##fileformat=VCFv4.1\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for var in variants:
                out_f.write(f"{var['chrom']}\t{var['pos']}\t"
                            f".\t{var['ref']}\t{var['alt']}\t.\t.\t.\n")
