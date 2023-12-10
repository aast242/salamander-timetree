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
import os
import itertools
import sys
from shutil import rmtree
from pprint import pprint

import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Salamander Taxfile\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("scraped_genes_fasta", help='File containing genes from gene_scraper.py', type=str)

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delieted on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()


if __name__ == '__main__':
    seq_records = list(SeqIO.parse(args.scraped_genes_fasta, "fasta"))

    tax_names = set(k.description.replace(k.id, "").replace(k.id.split(".")[-1], "").strip() for k in seq_records)
    tax_names = set(k.split("DESCRIPTION")[0].strip() for k in tax_names)

    # update sequence dictionary
    k = 0
    while k < len(seq_records):
        temp_sp = seq_records[k].description.replace(seq_records[k].id, "")
        temp_sp = temp_sp.replace(seq_records[k].id.split(".")[-1], "").strip()
        temp_sp = temp_sp.replace("DESCRIPTION", "").strip()
        temp_sp = temp_sp.replace("*", "").strip()
        temp_sp = temp_sp.replace(" ", "-")
        seq_records[k].id = "%s.%s" % (temp_sp, seq_records[k].id)
        seq_records[k].name = "%s.%s" % (temp_sp, seq_records[k].name)
        seq_records[k].description = "%s.%s" % (temp_sp, seq_records[k].description)
        k += 1

    # write taxonomy files
    SeqIO.write(seq_records, "%s_spp-reformat.fasta" % Path(args.scraped_genes_fasta).stem, "fasta")
