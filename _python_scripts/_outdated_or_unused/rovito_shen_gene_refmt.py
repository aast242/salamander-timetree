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

from defaults import ProgDefaults as dv
import utils

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Rovito et al. and Shen et al. Gene Table Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("ncbi_accs", help='A tsv file containing NCBI accessions for the genes used in'
                                      'Rovito et al. 2015, and Shen et al. 2013 & 2016', type=str)
parser.add_argument("--scraped_genes", default="", type=str,
                    help="Path to a tsv file that contains sequence records scraped from NCBI. Compares the IDs to the"
                         "IDs in the combined Rovito & Shen datafile.")
parser.add_argument("--gbk_file", action="store_true", help="Interprets the initial file as a gbk file")


args = parser.parse_args()


if __name__ == '__main__':
    if args.gbk_file:
        genbank_records = SeqIO.parse(args.ncbi_accs, "genbank")
        for record in genbank_records:
            print(record.annotations)
        exit()
    raw_acc_data = utils.parse_file_nohead_tolist(args.ncbi_accs)
    # use a set to resolve any redundancy present in the Shen data
    unique_accs = set()

    # some of the rovito accessions are paired and separated with a comma
    for p in raw_acc_data:
        for k in p:
            comma_split = k.split(",")
            if len(comma_split) > 1:
                for part in comma_split:
                    unique_accs.add(part.rstrip())
            else:
                unique_accs.add(k)

    # remove the placeholder characters for gaps: "" and "—"
    unique_accs = [k.strip() for k in sorted(list(unique_accs)) if k != "" and k != "—"]

    if args.scraped_genes != "":
        print("Length before removing duplicate IDs: %s" % len(unique_accs))
        scraped_ids = utils.parse_file_nohead_tolist(args.scraped_genes)
        scraped_ids = set(k[dv.ORDER_DICT["accession"]].split(".")[0] for k in scraped_ids)
        unique_accs = [k for k in unique_accs if k not in scraped_ids]
        print("Length after removing duplicate IDs: %s" % len(unique_accs))

    handle = Entrez.efetch(db="nucleotide", id=unique_accs, rettype="gb", retmode="gbwithparts")
    with open("rovito_shen_novel.gbk", "w") as f:
        for i in list(handle):
            f.write(str(i))
