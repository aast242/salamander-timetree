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
import itertools
import sys
from shutil import rmtree
from pprint import pprint

import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Remove redundancy in fasta\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("target_fasta", help='', type=str)

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delieted on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()

if __name__ == '__main__':
    target_seqs = list(SeqIO.parse(args.target_fasta, "fasta"))


    SeqIO.write(target_seqs, "%s_nr.fasta" % Path(args.target_fasta).stem, "fasta")
