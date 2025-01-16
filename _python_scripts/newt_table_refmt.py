#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez, SeqIO, SeqRecord, Seq
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

parser = argparse.ArgumentParser(description="Program: Pyron Data Merging\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("table", help='A file containing newt sequences', type=str)

newt_tt = {"E.": "Echinotriton", "T.": "Tylototriton", "P.": "Pleurodeles"}

args = parser.parse_args()

if __name__ == '__main__':
    newt_table = utils.parse_file_nohead_tolist(args.table)
    del newt_table[0]
    i = 0
    while i < len(newt_table):
        dash_count = newt_table[i][1].count("-")
        # this is a range of values that will need to be expanded
        if dash_count < 1:
            pass
        # these are individual values that should each get their own line
        else:
            genbank_eles = newt_table[i][1].split("-")
            genbank_eles = [[newt_table[i][0], k] for k in genbank_eles]
            newt_table.extend(genbank_eles)
            del newt_table[i]
            continue
        i += 1
    newt_table = [i.split("$$$") for i in list(set(["%s$$$%s" % (k[0], k[1]) for k in newt_table]))]
    newt_table.sort(key=lambda k: k[0])
    for i in newt_table:
        if " sp. nov." in i[0]:
            i[0] = i[0].replace(" sp. nov.", "")
        sp_parts = i[0].split(" ")
        sp_parts[0] = newt_tt[sp_parts[0]]
        i[0] = " ".join(sp_parts)

    with open("newt_table_reformatted.txt", 'w') as f:
        for line in newt_table:
            for ele in line:
                f.write(ele)
                f.write("\t")
            f.write("\n")
