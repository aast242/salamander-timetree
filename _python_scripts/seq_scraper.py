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
import os
from pathlib import Path
from shutil import rmtree
from collections import Counter
from pprint import pprint

import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Salamander Sequence Scraper\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("scraped_genes", help='File containing genes from gene_scraper.py', type=str)
parser.add_argument("mode", choices=["scrape", "split"], type=str,
                    help="which ")
parser.add_argument("--force", help='forces deleting files', action='store_true')

args = parser.parse_args()


def scrape_sequences(chunked_acc_list):
    try:
        if args.force:
            try:
                rmtree(dv.OUTDIR_PREFIX)
            except FileNotFoundError:
                pass
        os.mkdir(dv.OUTDIR_PREFIX)
    except FileExistsError:
        print("FATAL: Cannot create \'%s\' directory. Use \'--force\' to overwrite files" % dv.OUTDIR_PREFIX)
        exit()

    # iterate through each accession chunk and download the sequence for it
    print("There are %s accession chunks!" % len(chunked_acc_list))
    for acc in range(0, len(chunked_acc_list)):

        filename = "%s/chunk_%s.fa" % (dv.OUTDIR_PREFIX, acc)
        print("fetching and writing chunk_%s..." % acc)

        handle = Entrez.efetch(db='nucleotide', id=chunked_acc_list[acc], rettype="fasta", retmode="text")

        with open(filename, "w") as f:
            f.write(handle.read())

        handle.close()


def combine_fasta_chunks(chunked_acc_list):
    seq_records = []
    for i in range(0, len(chunked_acc_list)):
        opened_chunk = SeqIO.parse("%s/chunk_%s.fa" % (dv.OUTDIR_PREFIX, i), "fasta")
        seq_records.extend(list(opened_chunk))
    return SeqIO.to_dict(seq_records)
    # SeqIO.write(seq_records, "%s/%s" % (OUTDIR_PREFIX, COMB_CHUNK_NAME), "fasta")


def make_gene_dirs(accessions):
    uniq_genes = list(set(i[dv.ORDER_DICT["gene"]] for i in accessions))
    for i in uniq_genes:
        try:
            os.mkdir("%s/%s" % (dv.OUTDIR_PREFIX, i))
        except FileExistsError:
            pass


def subsequence_files(accessions, sequences):
    #make_gene_dirs(accessions)

    entries_by_acc = {}
    renamed_records = []
    for i in accessions:
        try:
            entries_by_acc[i[dv.ORDER_DICT["accession"]]].append(i)
        except KeyError:
            entries_by_acc.setdefault(i[dv.ORDER_DICT["accession"]], [i])

    for i in entries_by_acc.keys():
        active_seq = sequences[i]
        for k in entries_by_acc[i]:
            start, stop = sorted([int(j) for j in k[dv.ORDER_DICT["range"]].split(":")])

            # modify indices to fit with python
            start -= 1
            stop = stop

            temp_subseq = active_seq[start:stop]

            # modify name
            temp_subseq.id += ".%s" % k[dv.ORDER_DICT["gene"]]
            temp_subseq.description = "%s %s" % (k[dv.ORDER_DICT["species"]], k[dv.ORDER_DICT["gene"]])

            renamed_records.append(temp_subseq)
    SeqIO.write(renamed_records, "%s/%s.fasta" % (dv.OUTDIR_PREFIX, dv.COMB_FILE), "fasta")


if __name__ == '__main__':
    g_path = Path(args.scraped_genes).resolve()
    g_list = utils.parse_file_nohead_tolist(g_path)

    # get accessions
    acc_list = list(set(k[dv.ORDER_DICT["accession"]] for k in g_list))
    acc_list.sort()

    chunked_list = list(utils.chunks(acc_list, dv.CHUNK_SIZE))

    #split_fasta_chunks(chunked_list)
    if args.mode == "split":
        print(len(g_list))
        utils.check_gaps(g_list)
        print(len(g_list))
        subsequence_files(g_list, combine_fasta_chunks(chunked_list))
    elif args.mode == "scrape":
        scrape_sequences(chunked_list)
