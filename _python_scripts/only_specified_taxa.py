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

parser = argparse.ArgumentParser(description="Program: Only Specified Taxa\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("gene_fasta", help="", type=str)
parser.add_argument("exclude_file", help="", type=str)
args = parser.parse_args()


if __name__ == '__main__':
    file_stem = Path(args.gene_fasta).stem
    # make a directory to store outputs
    try:
        os.mkdir("Excluded_Fasta_Files")
    except FileExistsError:
        pass

    keep_names = [i[0] for i in utils.parse_file_nohead_tolist(args.exclude_file)]

    gene_seqs = SeqIO.to_dict(SeqIO.parse(args.gene_fasta, "fasta"))
    keep_dict = {}
    for ident in keep_names:
        try:
            keep_dict[ident] = gene_seqs[ident]
        except KeyError:
            pass

    gene_seqs = list(keep_dict.values())
    SeqIO.write(gene_seqs, f"Excluded_Fasta_Files/{file_stem}_SeqExclude.fasta", "fasta")
