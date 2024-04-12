#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez, SeqIO
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime
import os
import itertools
import sys
from shutil import rmtree
from pprint import pprint


parser = argparse.ArgumentParser(description="Program: Make supermatrix fig\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("captured_taxa", help='File containing captured taxa', type=str)
parser.add_argument("gene_dir", help='directory containing fasta files w/ genes', type=str)


# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delieted on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()


def parse_file_nohead_tolist(foi):
    parse_list = []
    cycle = 0
    with open(foi, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = line.rstrip().split('\t')
            parse_list.append(line)
            cycle += 1
    return parse_list


if __name__ == '__main__':
    taxon_names = [i[0] for i in parse_file_nohead_tolist(args.captured_taxa)]
    fasta_paths = ["%s/%s" % (Path(args.gene_dir), i) for i in os.listdir(args.gene_dir) if i.endswith(".fasta")]
    gene_names = [i.split("oneseq")[0].strip("_") for i in os.listdir(args.gene_dir) if i.endswith(".fasta")]
    fasta_dict = {}
    matrix_dict = {i: [] for i in gene_names}
    for i in range(0, len(fasta_paths)):
        fasta_dict[fasta_paths[i]] = gene_names[i]
    for i in fasta_paths:
        fasta_records = SeqIO.to_dict(SeqIO.parse(i, "fasta"))
        temp_entry = []
        for k in taxon_names:
            try:
                fasta_records[k]
                temp_entry.append(1)
            except KeyError:
                temp_entry.append(0)
        matrix_dict[fasta_dict[i]] = temp_entry

    final_matrix = []
    for i in matrix_dict.keys():
        final_matrix.append(matrix_dict[i])
    long_matrix = []
    for i in range(0, len(final_matrix)):
        for k in range(0, len(final_matrix[i])):
            long_matrix.append([gene_names[i], taxon_names[k], final_matrix[i][k]])

    with open("supermatrix_long.tsv" , 'w') as f:
        for x in long_matrix:
            for k in x:
                f.write("%s" % k)
                f.write("\t")
            f.write("\n")
