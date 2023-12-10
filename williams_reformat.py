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

parser = argparse.ArgumentParser(description="Program: Williams et al. Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("williams_nexus_dir", help='A *DIRECTORY* containing all of the williams nexus files', type=str)
# this file was retrieved from the dryad repository, subfolder 'MrBayes nexus Files'
# mitochondrial genes were removed (already sampled via other searches) and only the _wd files were kept
# for file pairs that had two. The _wd suffix seems to mean "with Dicamptodon". More sampling is good!
parser.add_argument("translation_table", type=str,
                    help='A taxonomy translation table from what is in the nexus files to latin names')
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")

args = parser.parse_args()

# translation table for homologous genes to Everson et al. 2021 determined via blast
EVERSON_GENES = {
    "cd163l1": "e24a12",
    "fmo3": "e24a6",
    "iqgap1": "e12g2",
    "lhx2": "e26g9",
    "sec22b": "e20c2",
    "trmt5": "e12c11"
}

if __name__ == '__main__':
    # set up a folder for exports
    nexus_dir_name = Path(args.williams_nexus_dir).stem
    export_dir_name = Path(f"./{nexus_dir_name}_SC-reformat")
    try:
        os.mkdir(export_dir_name)
    except FileExistsError:
        pass

    # make a translation table for species names
    dww_tt = utils.parse_file_nohead_todict(args.translation_table, 0)
    dww_tt = {i: dww_tt[i][1] for i in dww_tt.keys()}

    # get the names of all of the nexus files to modify
    nexus_files = [Path(f"{Path(args.williams_nexus_dir)}/{i}") for i in os.listdir(nexus_dir_name)
                   if i.endswith("nex")]

    # iterate through each nexus file
    for file_name in nexus_files:
        # get the gene name from the nexus file name
        gene_name = file_name.with_suffix("").stem.lower()
        # remove rubbish
        gene_name = gene_name.replace("ap_", "").replace("_wd", "")

        if gene_name in EVERSON_GENES.keys():
            print(gene_name)
            gene_name = EVERSON_GENES[gene_name]
            print(gene_name)
            print("----------")

        # open the nexus file
        gene_dictionary = SeqIO.to_dict(SeqIO.parse(file_name, "nexus"))
        fasta_ids = list(gene_dictionary.keys())
        for entry in fasta_ids:

            # If there are Dicamptodon sequences, get a different part of the name
            if entry.split("_")[0] == "D":
                sp_name = dww_tt[entry.split("_")[1]]
            else:
                sp_name = dww_tt[entry.split("_")[0]]

            # get the uid for the sequence (includes haplotype A/B)
            sp_uid = entry.split("_")[-1]

            # ungap the sequence
            gene_dictionary[entry].seq = gene_dictionary[entry].seq.ungap("-")
            # change undetermined characters from ? to N
            gene_dictionary[entry].seq = Seq.Seq(str(gene_dictionary[entry].seq).replace("?", "N"))

            # reformat the parts of the entries
            gene_dictionary[entry].id = f"{sp_uid}.{gene_name}"
            gene_dictionary[entry].name = gene_dictionary[entry].id
            gene_dictionary[entry].description = f"{sp_name} {gene_name}"
        gene_fasta = gene_dictionary.values()
        uniq_names = set()
        for entry in gene_fasta:
            uniq_names.add(" ".join(entry.description.split(" ")[:2]))
        if not args.quiet:
            pass
            #print(sorted(list(uniq_names)))
        SeqIO.write(gene_fasta, f"{export_dir_name}/{gene_name}_refmt.fasta", "fasta")
