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
parser.add_argument("out6dir", help='A *DIRECTORY* containing all of the blast outfiles', type=str)
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")
parser.add_argument("--hime", action="store_true", help="Parse things like the db is hime")
parser.add_argument("--split_char", default=".", type=str,
                    help="what character the gene name is split at")
parser.add_argument("--evalue", default=1e-10, type=float,
                    help="lowest e-value to consider hits for")

args = parser.parse_args()


def parse_out6(foi):
    parse_list = []
    cycle = 0
    with open(foi, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = line.rstrip().split('\t')
            for i in range(0, len(line)):
                line[i] = dv.BLASTN_TYPES[dv.BLASTN_DEFAULT_OUTFMT_SIX[i]](line[i])
            parse_list.append(line)
            cycle += 1
    return parse_list


if __name__ == '__main__':
    # get the names of all of the out6 files to parse
    data_files = [Path(f"{Path(args.out6dir)}/{i}") for i in os.listdir(Path(args.out6dir))
                  if i.endswith("out6")]
    # iterate through each nexus file
    translation_dictionary = {}
    for file_name in data_files:
        curr_file = parse_out6(file_name)
        curr_file = [k for k in curr_file if k[10] <= args.evalue]

        if len(curr_file) > 10:
            longest_aln = max([k[3] for k in curr_file])
            #if args.hime:
                #long_enough_ids = [k for k in curr_file if k[3] >= (longest_aln * 0.75)]
                #long_enough_ids = [k[1].split(args.split_char)[-1] for k in long_enough_ids if k[2] >= 80.00 and k[3] >= 150]
            #else:
            long_enough_ids = [k[1].split(args.split_char)[-1] for k in curr_file if k[3] >= (longest_aln * 0.75)]
            if len(long_enough_ids) == 0:
                continue

            most_common_id = Counter(long_enough_ids)
            if len(most_common_id) > 1:
                print(file_name)
                print(most_common_id)
            most_common_id = most_common_id.most_common(1)[0][0]
            translation_dictionary[str(Path(file_name).stem).replace("_refmt.fasta", "")] = most_common_id
            print("%s: %s" % (str(file_name.stem).split("_")[0], most_common_id))
    with open("temp_mergefile", 'w') as f:
        for locusname_2022 in translation_dictionary.keys():
            f.write(locusname_2022)
            f.write("\t")
            f.write(translation_dictionary[locusname_2022])
            f.write("\n")
