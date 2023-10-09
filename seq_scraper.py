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
from re import findall
import itertools
import sys
from shutil import rmtree
from collections import Counter
from pprint import pprint
import time

Entrez.email = "SETME"
Entrez.api_key = "SETME"

parser = argparse.ArgumentParser(description="Program: Salamander Sequence Scraper\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("scraped_genes", help='File containing genes from gene_scraper.py', type=str)
parser.add_argument("mode", choices=["scrape", "split"], type=str,
                    help="which ")
parser.add_argument("--force", help='forces deleting files', action='store_true')

args = parser.parse_args()

ORDER_DICT = {"species": 0, "accession": 1, "range": 2, "gene": 3, "date": 4, "length": 5}
OUTDIR_PREFIX = "scraped_sequences"
COMB_FILE = "all_chunks.fasta"
CHUNK_SIZE = 2000


def scrape_sequences(chunked_acc_list):
    try:
        if args.force:
            try:
                rmtree(OUTDIR_PREFIX)
            except FileNotFoundError:
                pass
        os.mkdir(OUTDIR_PREFIX)
    except FileExistsError:
        print("FATAL: Cannot create \'%s\' directory. Use \'--force\' to overwrite files" % OUTDIR_PREFIX)
        exit()

    # iterate through each accession chunk and download the sequence for it
    print("There are %s accession chunks!" % len(chunked_acc_list))
    for acc in range(0, len(chunked_acc_list)):

        filename = "%s/chunk_%s.fa" % (OUTDIR_PREFIX, acc)
        print("fetching and writing chunk_%s..." % acc)

        handle = Entrez.efetch(db='nucleotide', id=chunked_acc_list[acc], rettype="fasta", retmode="text")

        with open(filename, "w") as f:
            f.write(handle.read())

        handle.close()


def combine_fasta_chunks(chunked_acc_list):
    seq_records = []
    for i in range(0, len(chunked_acc_list)):
        opened_chunk = SeqIO.parse("%s/chunk_%s.fa" % (OUTDIR_PREFIX, i), "fasta")
        seq_records.extend(list(opened_chunk))
    return SeqIO.to_dict(seq_records)
    # SeqIO.write(seq_records, "%s/%s" % (OUTDIR_PREFIX, COMB_CHUNK_NAME), "fasta")


def make_gene_dirs(accessions):
    uniq_genes = list(set(i[ORDER_DICT["gene"]] for i in accessions))
    for i in uniq_genes:
        try:
            os.mkdir("%s/%s" % (OUTDIR_PREFIX, i))
        except FileExistsError:
            pass


def subsequence_files(accessions, sequences):
    #make_gene_dirs(accessions)

    entries_by_acc = {}
    renamed_records = []
    for i in accessions:
        try:
            entries_by_acc[i[ORDER_DICT["accession"]]].append(i)
        except KeyError:
            entries_by_acc.setdefault(i[ORDER_DICT["accession"]], [i])

    for i in entries_by_acc.keys():
        active_seq = sequences[i]
        for k in entries_by_acc[i]:
            start, stop = sorted([int(j) for j in k[ORDER_DICT["range"]].split(":")])

            # modify indices to fit with python
            start -= 1
            stop = stop

            temp_subseq = active_seq[start:stop]

            # modify name
            temp_subseq.id += ".%s" % k[ORDER_DICT["gene"]]
            temp_subseq.description = "%s %s" % (k[ORDER_DICT["species"]], k[ORDER_DICT["gene"]])

            renamed_records.append(temp_subseq)
    SeqIO.write(renamed_records, "%s/%s" % (OUTDIR_PREFIX, COMB_FILE), "fasta")


def parse_file_nohead(foi):
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
            parse_list.append(line)
            cycle += 1
    return parse_list


def check_gaps(gene_list):
    accessions = list(set(k[ORDER_DICT["accession"]] for k in gene_list))
    accessions.sort()

    ungapped_entries = []
    orig_gapped = []
    grouped_accessions = {i: [] for i in accessions}
    for gene in gene_list:
        grouped_accessions[gene[ORDER_DICT["accession"]]].append(gene)

    for i in grouped_accessions.keys():
        acc_count = Counter([j[ORDER_DICT["gene"]] for j in grouped_accessions[i]])
        max_num_genes = acc_count.most_common()[0][1]
        if max_num_genes > 1:
            gapped_genes = [j[0] for j in acc_count.most_common() if j[1] == max_num_genes]
            for gene in gapped_genes:
                more_than_two = [j for j in grouped_accessions[i] if j[ORDER_DICT["gene"]] == gene]

                range_nums = []
                for j in more_than_two:
                    for k in j[ORDER_DICT["range"]].split(":"):
                        range_nums.append(int(k))
                range_nums.sort()

                new_entry = more_than_two[0].copy()
                new_entry[ORDER_DICT["range"]] = "%s:%s" % (min(range_nums), max(range_nums))
                new_entry[ORDER_DICT["length"]] = "%s" % (max(range_nums) - min(range_nums))
                print(new_entry)
                ungapped_entries.append(new_entry)

                for j in more_than_two:
                    orig_gapped.append(j)

    for i in orig_gapped:
        gene_list.remove(i)
    for i in ungapped_entries:
        gene_list.append(i)

    gene_list.sort(key=lambda i: i[ORDER_DICT["species"]])


# from https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
# authored by Ned Batchelder
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


if __name__ == '__main__':
    g_path = Path(args.scraped_genes).resolve()
    g_list = parse_file_nohead(g_path)

    # get accessions
    acc_list = list(set(k[ORDER_DICT["accession"]] for k in g_list))
    acc_list.sort()

    chunked_list = list(chunks(acc_list, CHUNK_SIZE))

    #split_fasta_chunks(chunked_list)
    if args.mode == "split":
        print(len(g_list))
        check_gaps(g_list)
        print(len(g_list))
        subsequence_files(g_list, combine_fasta_chunks(chunked_list))
    elif args.mode == "scrape":
        scrape_sequences(chunked_list)
