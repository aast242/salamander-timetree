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

parser = argparse.ArgumentParser(description="Program: Rancilhac et al. Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("rancilhac_aln_dir", help='A *DIRECTORY* containing all of the Rancilhac fastas', type=str)
parser.add_argument("trans_table", help="A file containing the translation from Rancilhac gene names to main"
                                        "dataset names")
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")

# >KEBL01000001.1.prdm16 Caecilia tentaculata prdm16

species_to_change = {"Lyciasalamandra luschani": "Lyciasalamandra luschani",
                     "Lyciasalamandra basoglui": "Lyciasalamandra luschani",
                     "Lyciasalamandra finikensis": "Lyciasalamandra luschani",
                     "Lyciasalamandra flavi": "Lyciasalamandra flavimembris"
                     }
args = parser.parse_args()


def refmt_aln_file(aln_path, trans_table):
    gene_name = aln_path.with_suffix("").stem.lower().replace("-", "_")
    gene_name = trans_table[gene_name]
    aln_list = utils.parse_file_nohead_tolist(aln_path)
    for line in aln_list:
        if line[0].startswith(">"):
            # get the name of the species from fasta header lines
            sp_name = line[0].split("@")[0].split("_")[0].replace(">", "").strip()
            if sp_name in species_to_change.keys():
                sp_name = species_to_change[sp_name]
            # get the UID of the sequence from the fasta header lines and replace . with _
            sq_uid = line[0].split("@")[1].replace(".", "_").lower()

            # print(f"gene: {gene_name}\nsp: {sp_name}\nseq_uid: {sq_uid}\n\n\n")
            new_seqline = f">{sq_uid}.{gene_name} {sp_name} {gene_name}"
            line[0] = new_seqline
        else:
            # remove spaces and asterisks from alignments
            stripped_seq = line[0].replace(" ", "").replace("*", "")
            line[0] = stripped_seq
    return aln_list


if __name__ == '__main__':
    # set up the translation table
    tt = utils.parse_file_nohead_todict(args.trans_table, 0)
    tt = {i: tt[i][1] for i in tt.keys()}

    # set up a folder for exports
    aln_dir_name = Path(args.rancilhac_aln_dir).stem
    export_dir_name = Path(f"./{aln_dir_name}_sc-reformat")
    try:
        os.mkdir(export_dir_name)
    except FileExistsError:
        pass

    # get the names of all of the fasta files to modify
    newt_alns = [Path(f"{Path(args.rancilhac_aln_dir)}/{i}") for i in os.listdir(aln_dir_name)
                 if i.endswith("ali")]

    # find which fasta files are valid genes to include in the matrix
    newt_alns = [i for i in newt_alns if i.with_suffix("").stem.lower().replace("-", "_") in tt.keys()]

    # iterate through each valid fasta file
    for file_name in newt_alns:
        # get the gene name from the fasta file name
        gene_name = file_name.with_suffix("").stem.lower().replace("-", "_")
        gene_name = tt[gene_name]
        # open the fasta file
        gene_fasta = refmt_aln_file(file_name, tt)
        pprint([i[0] for i in gene_fasta if "Lyciasalamandra flavi" in i[0]])

        with open(f"{export_dir_name}/{gene_name}_refmt.fasta", "w") as f:
            for line in gene_fasta:
                f.write(line[0])
                f.write("\n")
