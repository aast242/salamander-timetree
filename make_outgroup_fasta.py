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
import os
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
parser.add_argument("outgroup_genes", type=str,
                    help='path to combined outgroup file made using \'outgroup_pull_genes.py\'')

args = parser.parse_args()


if __name__ == '__main__':
    outgroup_dict = utils.parse_file_head_todict(args.outgroup_genes, 0)
    ncbi_acc_dict = {}
    idx_dict = {k: outgroup_dict["header"][k] for k in range(0, len(outgroup_dict["header"]))}
    for sp in outgroup_dict.keys():
        if sp == "header":
            continue
        for element in range(1, len(outgroup_dict[sp])):
            ele_ref = outgroup_dict[sp][element]

            temp_ele = ele_ref.split(":")
            if len(temp_ele) == 1:
                continue
            ncbi_acc_dict.setdefault(temp_ele[0], [])

            temp_ele[1] = [int(k) for k in sorted(temp_ele[1].split("-"))]

            if temp_ele[2] == "plus":
                seq_pol = +1
            else:
                seq_pol = -1

            element_feature = SeqFeature.FeatureLocation(min(temp_ele[1])-1, max(temp_ele[1]), strand=seq_pol)
            ncbi_acc_dict[temp_ele[0]].append([element_feature, idx_dict[element], sp])

    output_dir = "outgroup_chroms"
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    goi_records = []
    loop_iter = 1
    for chrom in ncbi_acc_dict.keys():
        if not Path("%s/%s.fasta" % (output_dir, chrom)).exists():
            print("Fetching sequence for %s" % chrom)
            handle = Entrez.efetch(db="nucleotide", id=chrom, rettype="fasta", retmode="text")
            with open("%s/%s.fasta" % (output_dir, chrom), 'w') as f:
                for k in list(handle):
                    f.write(str(k))
        chrom_seq = list(SeqIO.parse("%s/%s.fasta" % (output_dir, chrom), "fasta"))[0].seq

        print("Getting GOIs from %s" % chrom)
        for goi in ncbi_acc_dict[chrom]:
            goi_records.append(SeqIO.SeqRecord(
                goi[0].extract(chrom_seq),
                id="%s.%s" % (chrom, goi[1]),
                name="%s.%s" % (chrom, goi[1]),
                description="%s %s" % (goi[2], goi[1])
                ))
        print("Done with %s! (%s/%s)\n" % (chrom, loop_iter, len(ncbi_acc_dict)))
        loop_iter += 1

    SeqIO.write(goi_records, "outgroup_gene_seqs.fasta", "fasta")
