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
import subprocess
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

parser = argparse.ArgumentParser(description="Program: Pull Gene from Genome\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("organism_genes", help='path to a tsv gene file from the genome of interest', type=str)
parser.add_argument("--add_dict", type=str, default="",
                    help="a file containing translations for gene names made using the \'--gbk_dir\' option"
                         " and !!CHECKED FOR INCONSISTENCIES!!")
parser.add_argument("--silent", help='do not print any messages', action='store_true')
# gene tsv can be downloaded at: https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000241765.5/
# get other sets by switching out the genbank assembly number

args = parser.parse_args()

if args.add_dict != "":
    utils.update_gene_dict_from_file(args.add_dict, dv.MITO_GOIS)


if __name__ == '__main__':
    ncbi_annotations = utils.parse_file_nohead_tolist(args.organism_genes)
    ann_order_dict = ncbi_annotations.pop(0)
    ann_order_dict = {ann_order_dict[k]: k for k in range(0, len(ann_order_dict))}
    gois = {k.lower(): [] for k in dv.MITO_GOIS.values()}

    for line in ncbi_annotations:
        line[ann_order_dict["Name"]] = line[ann_order_dict["Name"]].split(" isoform")[0]
        try:
            gn_in_goi = dv.MITO_GOIS[line[ann_order_dict["Name"]].lower()].lower()
        except KeyError:
            gn_in_goi = ""

        try:
            gs_in_goi = dv.MITO_GOIS[line[ann_order_dict["Symbol"]].lower()].lower()
        except KeyError:
            gs_in_goi = ""

        if gs_in_goi != "":
            gois[gs_in_goi].append(line)
        elif gn_in_goi != "":
            gois[gn_in_goi].append(line)

    for gene in gois.keys():
        if len(gois[gene]) > 1:
            match_count = []
            for potential in gois[gene]:
                temp_count = 0
                try:
                    # caecilians have some wonky low-quality genes
                    renamed_gn = potential[ann_order_dict["Name"]].replace("LOW QUALITY PROTEIN:", "").strip().lower()

                    if dv.MITO_GOIS[renamed_gn].lower() == gene:
                        temp_count += 1
                except KeyError:
                    pass
                try:
                    if dv.MITO_GOIS[potential[ann_order_dict["Symbol"]].lower()].lower() == gene:
                        temp_count += 1
                except KeyError:
                    pass
                match_count.append(temp_count)
            best_gois = [gois[gene][k] for k in range(0, len(match_count)) if match_count[k] == max(match_count)]
            best_gois.sort(key=lambda x: int(x[ann_order_dict["Protein length"]]), reverse=True)

            gois[gene] = [best_gois[0]]

        if len(gois[gene]) != 1:
            if not args.silent:
                print(len(gois[gene]))
                print(gene)
                print(gois[gene])
                print("--------")
            gois[gene] = ""
        else:
            gois[gene] = gois[gene][0]

            try:
                if gois[gene][ann_order_dict["Transcripts accession"]] != "":
                    new_id = gois[gene][ann_order_dict["Transcripts accession"]]
                else:
                    new_id = "%s:%s-%s:%s" % (gois[gene][ann_order_dict["Accession"]],
                                              gois[gene][ann_order_dict["Begin"]],
                                              gois[gene][ann_order_dict["End"]],
                                              gois[gene][ann_order_dict["Orientation"]])
            except IndexError:
                new_id = "%s:%s-%s:%s" % (gois[gene][ann_order_dict["Accession"]],
                                          gois[gene][ann_order_dict["Begin"]],
                                          gois[gene][ann_order_dict["End"]],
                                          gois[gene][ann_order_dict["Orientation"]])
            gois[gene] = new_id

    key_list = list(gois.keys())
    key_list.insert(0, "filename")
    value_list = list(gois.values())

    tsv_filename = str(Path(args.organism_genes).stem)
    tsv_filename = " ".join(tsv_filename.split("-", 1)[0].split("_"))
    value_list.insert(0, tsv_filename)

    with open("%s.genes" % Path(args.organism_genes).stem, 'w') as f:
        for line in [key_list, value_list]:
            for entry in line:
                f.write("%s" % entry)
                f.write("\t")
            f.write("\n")

