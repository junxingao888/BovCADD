#!/usr/bin/env python
# -*- coding: ASCII -*-
"""

:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 13.02.2018

- This file takes as input the different phyloFit models of each chunk and
computes for each substitution rate the weighted mean, it also reads in the
phylogenetic tree and computes the weighted mean of that.

- I need a dictionary with chromosome number as keys, values are
dictionaries with chunk numbers as keys and values are dictionaries with
path to model files, model file names, path to chunks, chunks file names,
chunk file size
- I need to know the total chunk size because in the end I
will have one model

:Edited by: Job van Schipstal
:Date: 4-10-2023
:Usage: see <script>.py --help
- Reformatted script to better follow PEP8 guidelines
- Moved from optparse to argparse
- Changed script to expect input model and chunks in subfolders for each chrom:
  chr[1]/[part].[mod/maf]
- Rewrote tree parsing with regex, now supports sequence names with numbers
  (e.g. previously mm39:0.6 would be parsed as mm=39 instead of mm39=0.6)
"""
import os
import re
import sys
from argparse import ArgumentParser

import numpy as np

TREE_WEIGHT_PATTERN = re.compile(r":([\d.]+)")

parser = ArgumentParser(description=__doc__)
parser.add_argument("-m", "--model-path",
                    help="Path to phyloFit model files",
                    type=str, required=True)
parser.add_argument("-c", "--chunk-path",
                    help="Path to phyloFit model files",
                    type=str, required=True)
parser.add_argument("-o", "--output",
                    help="Path to alignment chunks. This is necessary to get "
                         "the sizes of the files to compute the weighted "
                         "averages.",
                    type=str, required=True)


def path_check(path) -> str:
    """
    this function checks if input paths are non empty and end with an /
    it will return the input str object or a modified one
    :param path: str, path to be checked
    :return: str, path, always ending in /
    """
    # Error query and correction if the path is empty
    if path == '':
        sys.exit('Path to model files is not given.')

    # making sure that the path ensd with '/'
    if not path.endswith('/'):
        return path + '/'
    else:
        return path


def find_files(path, prefix="", suffix="") -> list:
    """
    This function identifies files in given directory with provided prefix and
    suffix It returns a file list containing all files in that directory,
    if empty stops program
    :param path: str, folder/path to look for files in
    :param prefix: optional str, prefix to filter files on
    :param suffix: optional str, suffix to filter files on
    :return: list of str, matched files, will exit if no files are found.
    """
    result_list = []
    for file in os.listdir(path):
        if file.startswith(prefix) and file.endswith(suffix):
            result_list.append(file)
    if len(result_list) == 0:
        sys.exit(
            'No files found in stated directory with given prefix and suffix '
            '\n Directory: %s \n Prefix: %s \n Suffix: %s' % (path, prefix,
                                                              suffix))
    else:
        return result_list


def dictionary_fill(path_models, path_chunks) -> (dict, int):
    """
    This function will create a dictionary with chromosome number as keys,
    values are dictionaries with chunk numbers as keys and values are
    dictionaries with path to model files, model file names, path to chunks,
    chunks file names, chunk file size
    :param path_models: str, path to model folder (.mod)
    :param path_chunks: str, path to chunk folder (.maf)
    :return: tuple of dict, described above and int, total size of chunk files.
    """
    # creating first layer of the dictionary by initializing dictionary with
    # chromosomes as keys
    outer_dict = {}
    folders = [name for name in os.listdir(path_models)
               if os.path.isdir(path_models + name) and name.startswith("chr")]
    if len(folders) == 0:
        sys.exit(f"No subfolders of chr[chr] found in {path_models}")
    for subfolder in folders:
        chrom = subfolder[3:]
        outer_dict[chrom] = {}

    for chrom in outer_dict.keys():
        # Creating second layer with chunks as keys and also adding the model
        # name and paths
        model_path_chr = f"{path_models}chr{chrom}/"
        chunk_path_chr = f"{path_chunks}chr{chrom}/"
        for model in find_files(f"{path_models}chr{chrom}/", suffix=".mod"):
            chunk_number = model.split('.')[0]
            outer_dict[chrom][chunk_number] = {}
            outer_dict[chrom][chunk_number]['model_path'] = model_path_chr
            outer_dict[chrom][chunk_number]['chunk_path'] = chunk_path_chr
            outer_dict[chrom][chunk_number]['model_name'] = model

        # adding the file size of the chunk files, file size in bytes
        for chunk in find_files(f"{path_chunks}chr{chrom}/", suffix=".maf"):
            chunk_number = chunk.split('.')[0]
            outer_dict[chrom][chunk_number]['chunk_name'] = chunk
            # checking size of chunk file in Bytes
            outer_dict[chrom][chunk_number]['chunk_size'] = os.stat(
                outer_dict[chrom][chunk_number]['chunk_path'] +
                outer_dict[chrom][chunk_number]['chunk_name']).st_size

    # adding the summed file sizes per chromosome to a new dictionary,
    # this is necessary to compute the weighted average per chromosome,
    # creating a new dictionary is helpful when I want to iterate over the
    # other dict
    total_size = 0
    for chrom in outer_dict.keys():
        for chunk_number in outer_dict[chrom].keys():
            total_size += outer_dict[chrom][chunk_number]['chunk_size']
    return outer_dict, total_size


def computing_weights(outer_dict, total_size) -> dict:
    """
    This function will go through each lower level of the chromosomes and
    compute the weight which each chunk has of the total chromosome size
    :param outer_dict: dict in which to add chunk weights from size
    :param total_size: int, total size of all chunks
    :return: dict, outer_dict with chunk_sizes added
    """
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            outer_dict[chro][chunk_number]['chunk_weight'] = \
                outer_dict[chro][chunk_number]['chunk_size'] / float(
                    total_size)

    return outer_dict


def parsing_models(outer_dict) -> (dict, str):
    """
    This function is parsing the model files and adds to the storage dict 3
    more keys ['TRAINING_LNL']: float, ['BACKGROUND']: [float,float,float,
    float], ['RATE_MAT'] numpy.array(4x4), ['TREE'] [with as many floats as
    there are in the trees, the order is always the same because the same
    topology was fitted] as a second value it returns the tree without numbers
    :param outer_dict: dict[chrom][chunk] to parse
    :return: input outer dict, str of tree without weights
    """
    tree_string = None
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            model_infile = open(outer_dict[chro][chunk_number]['model_path'] +
                                outer_dict[chro][chunk_number]['model_name'])
            for lines in model_infile:
                if lines.startswith('TRAINING_LNL:'):
                    outer_dict[chro][chunk_number]['TRAINING_LNL'] = float(
                        lines.strip().split(' ')[1])
                elif lines.startswith('BACKGROUND:'):
                    background = list(map(float, lines.strip().split(' ')[1:]))
                    outer_dict[chro][chunk_number]['BACKGROUND'] = background
                elif lines.startswith('RATE_MAT:'):
                    rate_matrix = np.zeros((4, 4))
                    # The next 4 lines are all values of the substitution
                    # matrix, I will go though them manually
                    for i in range(4):
                        rates = model_infile.readline().split()
                        print(i, rates)
                        rate_matrix[int(i), :] = list(map(float, rates))
                    outer_dict[chro][chunk_number]['RATE_MAT'] = rate_matrix
                elif lines.startswith('TREE:'):
                    # This splits the tree and remove all non numerical non
                    # dot characters, empty strings are removed from the
                    # list and then all the numbers are converted to float
                    tree_string = lines.split()[1]
                    tree = list(map(float, re.findall(TREE_WEIGHT_PATTERN,
                                                      tree_string)))
                    outer_dict[chro][chunk_number]['TREE'] = tree

    # the tree is always the same therefore I have to remove the numbers and
    # keep the tree only onces the last lines should contain tree
    tree = re.sub(TREE_WEIGHT_PATTERN, ":", tree_string)
    print("Plain tree: ")
    print(tree)

    return outer_dict, tree


def computing_weighted_means(outer_dict):
    """
    This function computes the weighted means and returns the model parameters
    :param outer_dict: dict[chrom][chunk] to calculate weighted means from
    :return: float averaged_training_lnl,
             list averaged_background,
             4x4 array averaged_rate_mat,
             list averaged_tree
    """
    # Computing average 'TRAINING_LNL'
    tmp_list = []
    tmp_weights = []
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            tmp_list.append(outer_dict[chro][chunk_number]['TRAINING_LNL'])
            # Here I am writing out the weights for this chromosome,
            # this has to be done only once therefore it will not appear
            # later on anymore
            tmp_weights.append(outer_dict[chro][chunk_number]['chunk_weight'])

    averaged_training_lnl = np.average(tmp_list, weights=tmp_weights)
    tmp_list = []

    # Computing average Background
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            # I am already multiplying the weights with the values
            tmp_list.append(
                np.array(outer_dict[chro][chunk_number]['BACKGROUND']) *
                outer_dict[chro][chunk_number]['chunk_weight'])

    # Here I am summing the results over columns and have the new average
    # background
    averaged_background = list(np.sum(np.array(tmp_list), axis=0))
    tmp_array = np.zeros((4, 4))

    # Computing average RATE_MAT
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            tmp_array += (outer_dict[chro][chunk_number]['RATE_MAT'] *
                          outer_dict[chro][chunk_number]['chunk_weight'])

    averaged_rate_mat = tmp_array
    tmp_list = []

    # Computing average tree
    trees_chr_1 = next(iter(outer_dict.values()))
    tree_length = len(next(iter(trees_chr_1.values()))['TREE'])
    for chro in outer_dict.keys():
        for chunk_number in outer_dict[chro].keys():
            unweighted_tree = np.array(outer_dict[chro][chunk_number]['TREE'])
            if len(unweighted_tree) != tree_length:
                print(f"Discarding incomplete tree from chr{chro} part "
                      f"{chunk_number}")
                continue
            tmp_list.append(np.array(unweighted_tree) *
                            outer_dict[chro][chunk_number]['chunk_weight'])

    averaged_tree = list(np.sum(np.array(tmp_list), axis=0))

    return averaged_training_lnl, averaged_background, \
           averaged_rate_mat, averaged_tree


# This function takes the averaged values and the first model file as an
# example to name the new model file. It does not return anything
def print_model(averaged_training_lnl, averaged_background, averaged_rate_mat,
                averaged_tree, tree, output) -> None:
    """
    This function takes the averaged values and writes them in a .mod file.
    :param averaged_training_lnl: float, training lnl
    :param averaged_background: list, background
    :param averaged_rate_mat: 4x4 array of float
    :param averaged_tree: list of average distance value per node in tree
    :param tree: str, tree without distance values
    :param output: str, output file with path
    :return: None, printed to output file
    """
    output = open(output, 'w')

    output.write('ALPHABET: A C G T \n')
    output.write('ORDER: 0\n')
    output.write('SUBST_MOD: REV\n')
    output.write('TRAINING_LNL: %s\n' % averaged_training_lnl)
    output.write('BACKGROUND: %s \n' % (
        re.sub("[^0-9. -]", "", str(averaged_background))))
    output.write('RATE_MAT:\n')
    output.write(
        '  %s \n' % (re.sub("[^0-9. -]", "", str(averaged_rate_mat[0, :]))))
    output.write(
        '  %s \n' % (re.sub("[^0-9. -]", "", str(averaged_rate_mat[1, :]))))
    output.write(
        '  %s \n' % (re.sub("[^0-9. -]", "", str(averaged_rate_mat[2, :]))))
    output.write(
        '  %s \n' % (re.sub("[^0-9. -]", "", str(averaged_rate_mat[3, :]))))
    # to do the tree and add the averaged numbers I have to iterate over
    # each character and everytime when an : appears increase a counter
    # after inserting the next value from my list of averaed tree values
    i = 0
    output_tree = ''
    for char in tree:
        if char == ':':
            output_tree += char + str(averaged_tree[i])
            i += 1
        else:
            output_tree += char

    output.write('TREE: %s \n' % output_tree)

    output.close()


if __name__ == '__main__':
    args = parser.parse_args()
    # Path error check
    args.model_path = path_check(args.model_path)
    args.chunk_path = path_check(args.chunk_path)

    # function returns a dictionary of dictionaries of dictionaries. The lowest
    # level contains ['chro']['chunk_number']['chunk_size', 'model_name',
    # 'model_path', 'chunk_name', 'chunk_path']
    data_dict, total_sizes = dictionary_fill(args.model_path,
                                             args.chunk_path)

    # computing the weight of each chunk
    data_dict = computing_weights(data_dict, total_sizes)

    # parsing the model files
    data_dict, tree_without_lengths = parsing_models(data_dict)

    # computing the weighted averages
    TRAINING_LNL_averaged, BACKGROUND_averaged, RATE_MAT_averaged, \
    TREE_averaged = computing_weighted_means(data_dict)

    # printing new model file
    print_model(TRAINING_LNL_averaged, BACKGROUND_averaged, RATE_MAT_averaged,
                TREE_averaged, tree_without_lengths, args.output)
