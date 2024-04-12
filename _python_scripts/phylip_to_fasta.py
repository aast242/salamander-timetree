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

args = parser.parse_args()

if __name__ == '__main__':
    curr_file = utils.parse_phylip(args.phylip_file)
    SeqIO.write(curr_file, "%s" % Path(args.phylip_file).stem, "fasta")
