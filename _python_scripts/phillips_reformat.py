#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez, SeqIO, SeqRecord
import os
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime
import itertools
import sys
from shutil import rmtree
from pprint import pprint

from defaults import ProgDefaults as dv
import utils

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Phillips et al. Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("phil_aln_dir", help='A *DIRECTORY* containing all of the Phillips fastas', type=str)
parser.add_argument("--trans_table", help="A file containing the translation from Phillips gene names to main"
                                          "dataset names")
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")


args = parser.parse_args()


def refmt_aln_file(aln_path, trans_table):
    gene_name = aln_path.with_suffix("").stem.lower().replace("-", "_")
    if trans_table != {}:
        gene_name = trans_table[gene_name]

    aln_list = utils.parse_file_nohead_tolist(aln_path)
    for line in aln_list:
        if line[0].startswith(">"):
            # get the name of the species from fasta header lines
            #sp_name = " ".join(line[0].split("_")[:2]).replace(">", "")
            sp_name = " ".join(line[0].split("_")[:2]).replace(">", "")

            # get the UID of the sequence from the fasta header lines
            if line[0].count("_") > 1:
                sq_uid = "".join([line[0].split("_")[-1].lower(), "_",
                                  sp_name.split(" ")[0][0], sp_name.split(" ")[1][:3]]).lower()
            else:
                sq_uid = "".join(["000_", sp_name.split(" ")[0][0], sp_name.split(" ")[1][:3]]).lower()

            # print(f"gene: {gene_name}\nsp: {sp_name}\nseq_uid: {sq_uid}\n\n\n")
            new_seqline = f">{sq_uid}.{gene_name} {sp_name} {gene_name}"
            line[0] = new_seqline
        else:
            # remove spaces and asterisks from alignments
            stripped_seq = line[0].replace("-", "").replace("N", "")
            line[0] = stripped_seq
    aln_list = [i[0] for i in aln_list]
    # remove probe sequences
    i = 0
    while i < len(aln_list):
        if "probe" in aln_list[i].lower():
            del aln_list[i+1]
            del aln_list[i]
            i -= 1
        i += 1

    return aln_list


if __name__ == '__main__':
    # set up the translation table
    tt = {}
    if args.trans_table is not None:
        tt = utils.parse_file_nohead_todict(args.trans_table, 0)
        tt = {i: tt[i][1] for i in tt.keys()}

    # set up a folder for exports
    aln_dir_name = Path(args.phil_aln_dir).stem
    export_dir_name = Path(f"./{aln_dir_name}_sc-reformat")
    try:
        os.mkdir(export_dir_name)
    except FileExistsError:
        pass

    # get the names of all of the fasta files to modify
    eurycea_alns = [Path(f"{Path(args.phil_aln_dir)}/{i}") for i in os.listdir(aln_dir_name)
                    if i.endswith(".fas")]

    # find which fasta files are valid genes to include in the matrix
    if args.trans_table is not None:
        eurycea_alns = [i for i in eurycea_alns if i.with_suffix("").stem.lower().replace("-", "_") in tt.keys()]

    # iterate through each valid fasta file
    for file_name in eurycea_alns:
        # get the gene name from the fasta file name
        gene_name = file_name.with_suffix("").stem.lower().replace("-", "_")
        if args.trans_table is not None:
            gene_name = tt[gene_name]
            # open the fasta file
        gene_fasta = refmt_aln_file(file_name, tt)

        with open(f"{export_dir_name}/{gene_name}_refmt.fasta", "w") as f:
            for line in gene_fasta:
                f.write(line)
                f.write("\n")
