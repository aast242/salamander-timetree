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
import itertools
import os
import sys
import shutil
from pprint import pprint
from collections import Counter
import time

from defaults import ProgDefaults as dv
import utils

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Salamander Gene Scraper\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("gene_file", help='File containing genes to search for', type=str)
parser.add_argument("gene_name")

args = parser.parse_args()


if __name__ == '__main__':
    seq_records = SeqIO.parse(args.gene_file, "fasta")
    focal_records = [i for i in seq_records if i.id.endswith(args.gene_name)]
    SeqIO.write(focal_records, "%s-%s.fasta" % (args.gene_name, Path(args.gene_file).stem), "fasta")
