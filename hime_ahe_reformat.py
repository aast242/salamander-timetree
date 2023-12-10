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

parser = argparse.ArgumentParser(description="Program: Hime et al. AHE Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("ahe_gbff", help='A GenBank formatted file containing the sequences'
                                     ' from Hime et al. (2021)', type=str)
# this file can be downloaded at the link below:
# https://www.ncbi.nlm.nih.gov/Traces/wgs/?display=download
parser.add_argument("focal_taxa", help='File containing the species to pull from the Hime et al. data', type=str)
parser.add_argument("--add_dict", type=str, default="",
                    help="a file containing translations for gene names made using the \'--gbk_dir\' option"
                         " and !!CHECKED FOR INCONSISTENCIES!!")

args = parser.parse_args()

if args.add_dict != "":
    utils.update_gene_dict_from_file(args.add_dict, dv.MITO_GOIS)

# make gene description dict
GENE_DESCRIPTION = {gn: [] for gn in sorted(list(set(dv.MITO_GOIS.values())))}
for long_gn in dv.MITO_GOIS.keys():
    GENE_DESCRIPTION[dv.MITO_GOIS[long_gn]].append(long_gn)
for gn in GENE_DESCRIPTION.keys():
    GENE_DESCRIPTION[gn] = sorted(GENE_DESCRIPTION[gn], key=len, reverse=True)[0]

# harmonize nomenclature across datasources
for k in dv.HIME_TO_SYMBOL.keys():
    try:
        pre_change = dv.HIME_TO_SYMBOL[k]
        dv.HIME_TO_SYMBOL[k] = dv.MITO_GOIS[dv.HIME_TO_SYMBOL[k]]
        if pre_change != dv.HIME_TO_SYMBOL[k]:
            print(pre_change)
            print(dv.HIME_TO_SYMBOL[k])
            print("--------------")
    except KeyError:
        print("%s wasn't a key!" % k)
exit()
# dv.HIME_TO_SYMBOL = {k: "HL%s_%s" % (k, dv.HIME_TO_SYMBOL[k]) for k in dv.HIME_TO_SYMBOL.keys()}
SYMBOL_TO_HIME = {dv.HIME_TO_SYMBOL[k]: k for k in dv.HIME_TO_SYMBOL.keys()}


if __name__ == '__main__':
    raw_hime_data = SeqIO.parse(args.ahe_gbff, "genbank")

    focal_species = [k[0] for k in utils.parse_file_nohead_tolist(args.focal_taxa)]

    """
    salamander_spp = set()
    for record in raw_hime_data:
        if "Caudata" in record.annotations["taxonomy"]:
            salamander_spp.add(record.annotations["organism"])
    salamander_spp = sorted(list(salamander_spp))
    with open("hime_salamander_spp.txt", 'w') as f:
        for sp in salamander_spp:
            f.write(sp)
            f.write("\n")
    """

    hime_rois = []
    seq_dict = {}
    for record in raw_hime_data:
        seq_name = record.description.replace("TLS: %s " % record.annotations["organism"], "")
        seq_name = " ".join(seq_name.split(";")[0].split(" ")[1:])
        seq_name = int(seq_name.replace("anchored hybrid enrichment locus ", ""))
        seq_name = dv.HIME_TO_SYMBOL[seq_name]
        seq_id = "%s.%s" % (record.id, seq_name)

        """
        seq_long_desc = "label not included in dictionary"
        try:
            seq_long_desc = GENE_DESCRIPTION[seq_name]
        except KeyError:
            pass
        """

        temp_record = SeqRecord.SeqRecord(record.seq,
                                          description="%s %s" %
                                                      (record.annotations["organism"],
                                                       seq_name),
                                          name=seq_id,
                                          id=seq_id)
        seq_dict.setdefault(seq_name, [])
        seq_dict[seq_name].append(temp_record)
        if record.annotations["organism"] in focal_species:
            hime_rois.append(temp_record)

    SeqIO.write(hime_rois, "hime_ahe_reformatted.fasta", "fasta")
    exit()
    reference_seq_dir = "hime_reference_seqs"
    try:
        os.mkdir("hime_reference_seqs")
    except FileExistsError:
        pass
    for gene in seq_dict.keys():
        SeqIO.write(seq_dict[gene], "%s/hime-ref-%s.fasta" % (reference_seq_dir, gene), "fasta")
