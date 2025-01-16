#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

# base python imports #
import argparse
import os
import itertools
import sys
from shutil import rmtree
from pprint import pprint
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime

# external package imports #
from ete3 import PhyloTree
from Bio import Entrez, SeqIO, SeqRecord

# local imports #
from defaults import ProgDefaults as dv
import utils

parser = argparse.ArgumentParser(description="Program: phylip to fasta\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s phylip_file [options]')
parser.add_argument("phylip_file", help="", type=str)
parser.add_argument("--keep_ext", help="", action="store_true")

args = parser.parse_args()

if __name__ == '__main__':
    with open(args.phylip_file) as f:
        first_line = f.readline()
        if first_line.startswith(">"):
            print("FATAL: %s is not a phylip file!" % Path(args.phylip_file).stem)
            exit()
    curr_file = utils.parse_phylip(args.phylip_file)
    if args.keep_ext:
        fn = Path(args.phylip_file).name
    else:
        fn = Path(args.phylip_file).stem

    SeqIO.write(curr_file, "%s" % fn, "fasta")
