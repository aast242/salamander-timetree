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

parser = argparse.ArgumentParser(description="Program: Supplementary Table Data Gatherer\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("comb_gbk", type=str,
                    help='A .gbk file containing all the sequences scraped from genbank')
parser.add_argument("gois", type=str, help="A file containing the genbank IDs of the gois in the gbk file",
                    default="")
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")

# >KEBL01000001.1.prdm16 Caecilia tentaculata prdm16

args = parser.parse_args()

filename_tt = {
    "everson_2021_markers_refmt": "everson",
    "hime_ahe_reformatted": "hime",
    "hime-rovito-shen-scraped_genes_GenBank_rename": "genbank",
    "outgroup_1e-10": "outgroup",
    "pyronCombined_AHE_refmt.rename": "pyron",
    "rancilhac_2021_refmt": "rancilhac",
    "williams_2013_markers_refmt": "williams"
}


def parse_gbk_entry(gbk_record, genes_of_interest):
    # TODO: based on these analyses, there are a TON of "misc features" that have 12s and 16s fragments
    # maybe eventually try to pull those out? probably not worth it, but could be interesting
    # print(gbk_record.features)
    feature_list = []
    for feature in gbk_record.features:
        try:
            product_identity = feature.qualifiers["product"][0].lower()
            try:
                # Make sure this is a gene we actually care about!
                feature_gene = genes_of_interest[product_identity]

                feature_acc = gbk_record.id
                feature_organism = gbk_record.annotations["organism"]
                feature_range = feature.location
                feature_range = "%s:%s" % (int(feature_range.start), int(feature_range.end))
                feature_date = gbk_record.annotations["date"]
                feature_seq = feature.extract(gbk_record.seq)
                feature_length = len(feature_seq)

                feature_list.append([feature_organism, feature_acc, feature_range, feature_gene,
                                     feature_date, feature_length])

            except KeyError:
                pass
        except KeyError:
            pass
    return feature_list


if __name__ == '__main__':
    # parse the gois file
    parsed_goi = [i[0] for i in utils.parse_file_nohead_tolist(args.gois)]

    # parse the gbk file
    parsed_gbk = list(SeqIO.parse(args.comb_gbk, "genbank"))

    # only get the gois that we care about
    gbk_dict = {}
    info_dict = {}
    for record in parsed_gbk:
        if record.name in parsed_goi:
            gbk_dict.setdefault(record.name, record)
            curr_spp = gbk_dict[record.name].annotations["organism"]
            curr_vou = ""
            for ele in gbk_dict[record.name].features:
                if ele.type == "source":
                    try:
                        curr_vou = ele.qualifiers["specimen_voucher"][0]
                    except KeyError:
                        try:
                            curr_vou = ele.qualifiers["isolate"][0]
                        except KeyError:
                            pass
            info_dict[record.name] = [record.name, curr_spp, curr_vou]

    with open("s3_info.tsv", "w") as f:
        for line in list(info_dict.values()):
            for ele in line:
                f.write(ele)
                f.write("\t")
            f.write("\n")
