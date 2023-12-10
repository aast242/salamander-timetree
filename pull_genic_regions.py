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

parser = argparse.ArgumentParser(description="Program: Pull Genic Regions from Genome\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("gene_tsv", help="", type=str)
parser.add_argument("subj_fasta", help="")
parser.add_argument("--mito_acc", default="", type=str,
                    help="")
parser.add_argument("--wiggle", default=1000)
args = parser.parse_args()


if __name__ == '__main__':
    gene_regions = utils.parse_file_nohead_tolist(args.gene_tsv)
    gene_regions = [x for x in gene_regions if "pseudogene" not in x[8]]

    # 'Accession', 'Begin', 'End', 'Chromosome', 'Orientation', 'Name', 'Symbol',
    # 'Gene ID', 'Gene Type', 'Transcripts accession', 'Protein accession', 'Protein length'
    refmt_genes = []
    for x in range(1, len(gene_regions)):
        line = gene_regions[x]
        chrom_id = line[0]
        genic_range = sorted([int(line[1]), int(line[2])])
        genic_range = [genic_range[0] - args.wiggle, genic_range[1] + args.wiggle]
        if genic_range[0] < 0:
            genic_range[0] = 0

        orientation = line[4]
        if orientation == "minus":
            orientation = -1
        else:
            orientation = +1

        gene_name = line[6]
        uniq_gene_id = "%s_%s_%s" % (chrom_id, gene_name, x)

        refmt_genes.append([uniq_gene_id,
                            chrom_id,
                            SeqFeature.FeatureLocation(genic_range[0], genic_range[1], strand=orientation)])

    print("Indexing subject fasta. This might take a while...")
    # subject_idx = SeqIO.index(args.subj_fasta, "fasta")
    subject_idx = SeqIO.to_dict(SeqIO.parse(args.subj_fasta, "fasta"))
    genic_seqs = []
    for line in refmt_genes:
        temp_genic_seq = line[2].extract(subject_idx[line[1]]).seq
        temp_id = line[0]
        genic_seqs.append(SeqIO.SeqRecord(
            temp_genic_seq,
            id=temp_id,
            name=temp_id,
            description=""))
    try:
        mito_seq = subject_idx[args.mito_acc]
        genic_seqs.append(mito_seq)
        print("Mitochondrial sequence successfully added!")
    except KeyError:
        if args.mito_acc == "":
            print("Mitochondria was not added!")
        else:
            print("Mitochondrial accession not present in provided genome...")
            print("Mitochondria was not added!")

    SeqIO.write(genic_seqs, "%s-genic_regions.fasta" % Path(args.gene_tsv).stem, "fasta")
