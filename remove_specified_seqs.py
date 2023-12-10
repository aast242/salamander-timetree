#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez, SeqIO, SeqFeature
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime
from shutil import rmtree
from pprint import pprint
import os
import subprocess
import itertools
import sys


import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Remove Specified Sequences\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("gene_fasta", help="", type=str)
parser.add_argument("exclude_file", help="", type=str)
parser.add_argument("--rm_outgroups", help="ignore exlude file and remove records containing outgroup names",
                    action="store_true")
args = parser.parse_args()


if __name__ == '__main__':
    file_stem = Path(args.gene_fasta).stem
    # make a directory to store outputs
    try:
        os.mkdir("Excluded_Fasta_Files")
    except FileExistsError:
        pass

    if args.rm_outgroups:
        exclude_names = ["Latimeria chalumnae", "Homo sapiens", "Xenopus tropicalis", "Anolis carolinensis",
                         "Chrysemys picta", "Gallus gallus", "Geotrypetes seraphini", "Microcaecilia unicolor",
                         "Rhinatrema bivittatum"]
    else:
        exclude_names = [i[0] for i in utils.parse_file_nohead_tolist(args.exclude_file)]

    gene_seqs = SeqIO.to_dict(SeqIO.parse(args.gene_fasta, "fasta"))
    for ident in exclude_names:
        try:
            del gene_seqs[ident]
            print(f"deleted {ident}")
        except KeyError:
            pass

    gene_seqs = list(gene_seqs.values())
    SeqIO.write(gene_seqs, f"Excluded_Fasta_Files/{file_stem}_SeqExclude.fasta", "fasta")
