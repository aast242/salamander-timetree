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

parser = argparse.ArgumentParser(description="Program: Hime et al. AHE Refresher\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("stale_fasta", help='A concatenated fasta file containing similarity clustering-trimmed '
                                        'sequences from SuperCRUNCH', type=str)
parser.add_argument("hime_reformatted", help='The original \'hime_ahe_reformatted.fasta\' that entered the pipeline',
                    type=str)

args = parser.parse_args()

if __name__ == '__main__':
    stale_seqs = list(SeqIO.parse(args.stale_fasta, "fasta"))
    fresh_seqs = SeqIO.to_dict(SeqIO.parse(args.hime_reformatted, "fasta"))
    for entry in stale_seqs:
        if entry.id.startswith("KEBL"):
            temp_fresh = fresh_seqs[entry.id]
            entry.seq = temp_fresh.seq

    SeqIO.write(stale_seqs, "%s_refreshed.fasta" % Path(args.stale_fasta).stem, "fasta")
