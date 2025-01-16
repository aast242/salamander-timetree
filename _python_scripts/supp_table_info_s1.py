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

parser = argparse.ArgumentParser(description="Program: Supplementary Table Data Gatherer\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("pre_filter_num_seqs", type=str,
                    help='A *DIRECTORY* containing pre-filtering by number of seqs fastas')
parser.add_argument("--starting_seq_dir", type=str, help="A *DIRECTORY* containing the starting sequences for SC",
                    default="")
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")

# >KEBL01000001.1.prdm16 Caecilia tentaculata prdm16

args = parser.parse_args()

filename_tt = {
    "everson_2021_markers_refmt": "everson",
    "hime_ahe_reformatted": "hime",
    "hime-rovito-shen-scraped_genes_GenBank_rename": "genbank",
    "outgroup_1e-10": "outgroup",
    "pyronCombined_AHE_refmt.rename": "pyron",
    "rancilhac_2021_refmt": "rancilhac",
    "williams_2013_markers_refmt": "williams"
}

if __name__ == '__main__':
    starting_seq_tt = {}
    if args.starting_seq_dir != "":
        for i in ["%s/%s" % (Path(args.starting_seq_dir), i) for i in list(os.listdir(args.starting_seq_dir))]:

            starting_seq_tt[filename_tt[Path(i).stem]] = sorted(list({record.id: record for record
                                                                      in SeqIO.parse(i, "fasta")}.keys()))
    # set up the translation table
    rename_files = ["%s/%s" % (Path(args.pre_filter_num_seqs), i)
                    for i in list(os.listdir(args.pre_filter_num_seqs))]

    rename_dict = {}
    for i in rename_files:
        gn = Path(i).stem.split("_oneseq")[0]
        spp = sorted(list(SeqIO.to_dict(SeqIO.parse(i, "fasta")).keys()))
        rename_dict[gn] = spp
    final_info_list = []
    for gene in rename_dict.keys():
        for uid in rename_dict[gene]:
            curr_dataset = ""
            for starter_file in starting_seq_tt.keys():
                if uid in starting_seq_tt[starter_file]:
                    curr_dataset = starter_file
                    break
            final_info_list.append([gene, uid, curr_dataset])

    #pprint(final_info_list[:100])
    with open("s1_info.tsv", "w") as f:
        for line in final_info_list:
            for ele in line:
                f.write(ele)
                f.write("\t")
            f.write("\n")
