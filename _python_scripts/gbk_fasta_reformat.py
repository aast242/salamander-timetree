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

parser = argparse.ArgumentParser(description="Program: GenBank Fasta Reformatter\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("gene_fasta", help='fasta file containing all GenBank scraped sequences', type=str)
parser.add_argument("translation_table", help='a translation table with \"accession voucher species\"')
parser.add_argument("--exclude", help="a file containing accessions to exclude", default="")

args = parser.parse_args()

raise_spp = {"Lissotriton vulgaris graecus": "Lissotriton graecus",
             "Lissotriton vulgaris lantzi": "Lissotriton lantzi",
             "Lissotriton vulgaris schmidtlerorum": "Lissotriton schmidtleri",
             "Neurergus strauchii barani": "Neurergus barani",
             "Taricha torosa sierrae": "Taricha sierrae",
             "Hydromantes ambrosii": "Speleomantes ambrosii",
             "Triton torosa": "Taricha torosa",
             "Iranodon persicus": "Paradactylodon persicus",
             "Hypselotriton jiaoren": "Cynops jiaoren",
             "Hypselotriton maguae": "Cynops maguae",
             "Hypselotriton yunnanensis": "Cynops yunnanensis",
             "Nototriton matama": "Nototriton picadoi",
             "Paramesotriton ermizhaoi": "Paramesotriton labiatus",
             "Paramesotriton guanxiensis": "Paramesotriton guangxiensis",
             "Plethodon longicrus": "Plethodon yonahlossee",
             "Plethodon oconaluftee": "Plethodon teyahalee"
             }


if __name__ == '__main__':
    seq_records = list(SeqIO.parse(args.gene_fasta, "fasta"))
    tt = utils.parse_file_nohead_tolist(args.translation_table)
    exclude_list = []
    if args.exclude != "":
        exclude_list = [i[0] for i in utils.parse_file_nohead_tolist(args.exclude)]
    tt = {i[0]: i[2] for i in tt}
    for i in seq_records:
        focal_sp = i.description.replace(i.id, "").replace(i.id.split(".")[-1], "").strip()

        # elevate some subspecies to species as defined in the dictionary in the start of the file
        if focal_sp in raise_spp.keys():
            i.description = i.description.replace(focal_sp, raise_spp[focal_sp])

        if i.id.split(".")[0] in tt.keys():
            #if focal_sp != tt[i.id.split(".")[0]]:
            #    print(i.description)
            i.description = i.description.replace(focal_sp, tt[i.id.split(".")[0]])
            #if focal_sp != tt[i.id.split(".")[0]]:
            #    print(i.description)
            #    print("~~~~~~~~~~~~~~~~~~")

    focal_records = [i for i in seq_records if i.id.split(".")[0] not in exclude_list]
    excluded_records = [i for i in seq_records if i.id.split(".")[0] in exclude_list]

    SeqIO.write(focal_records, "%s_GenBank_rename.fasta" % Path(args.gene_fasta).stem, "fasta")

