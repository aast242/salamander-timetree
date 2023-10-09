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
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime
import itertools
import sys
from shutil import rmtree
from pprint import pprint

Entrez.email = "SETME"
Entrez.api_key = "SETME"

parser = argparse.ArgumentParser(description="Program: Salamander Gene Chooser\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("scraped_genes", help='File containing genes from gene_scraper.py', type=str)

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delieted on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()


ORDER_DICT = {"species": 0, "accession": 1, "range": 2, "gene": 3, "date": 4, "length": 5}
OUTDIR_PREFIX = "scraped_sequences"
COMB_FILE = "all_chunks.fasta"
EXCLUDE_TERMS = ["unisexual", "cf.", " x "]
SEQUENCE_DICTIONARY = SeqIO.to_dict(SeqIO.parse("%s/%s" % (OUTDIR_PREFIX, COMB_FILE), "fasta"))
FOCAL_SEQFILE_NAME = "focal_scraped_sequences.fasta"
SPP_FILE_NAME = "captured_caudate_taxa"


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
                # print(new_entry)
                ungapped_entries.append(new_entry)

                for j in more_than_two:
                    orig_gapped.append(j)

    for i in orig_gapped:
        gene_list.remove(i)
    for i in ungapped_entries:
        gene_list.append(i)

    gene_list.sort(key=lambda i: i[ORDER_DICT["species"]])


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


def get_uniq_element(parsed_file, element):
    return sorted(list(set(i[ORDER_DICT[element]] for i in parsed_file)))


def clean_species(spp_list):
    removed_spp = []
    i = 0
    while i < len(spp_list):
        #print([j in spp_list[i] for j in EXCLUDE_TERMS])
        if any([j in spp_list[i] for j in EXCLUDE_TERMS]):
            removed_spp.append(spp_list.pop(i))
        else:
            i += 1
    return removed_spp


def get_weird_spp(spp_list):
    novel_spp = []
    subspp = []
    for entry in spp_list:
        if len(entry.split(" ")) > 2:
            if "." in entry:
                novel_spp.append(entry)
            else:
                subspp.append(entry)
    """
    for entry in novel_spp:
        print(entry)
    print("------------------------------------")
    for entry in subspp:
        print(entry)
    """

    return novel_spp, subspp


def make_gene_dict(spp_list, gene_names, gois):
    final_dict = {i: {j: [] for j in gene_names} for i in spp_list}
    for entry in gois:
        try:
            final_dict[entry[ORDER_DICT["species"]]][entry[ORDER_DICT["gene"]]].append(entry)
        except KeyError:
            pass
    return final_dict


def filter_gene_dict(gene_dict):
    same_ct = 0
    diff_ct = 0

    for sp in gene_dict.keys():
        for gene in gene_dict[sp].keys():
            if len(gene_dict[sp][gene]) > 1:
                refseq_entries = ["_" in i[ORDER_DICT["accession"]] for i in gene_dict[sp][gene]]
                refseq_indices = [i for i in range(0, len(refseq_entries)) if refseq_entries[i] is True]

                # if there are refseq entries, those should be the representative sequence
                if any(refseq_entries):
                    # if there is more than one refseq entry (shouldn't happen, but added a case for safety)
                    # sort by length and choose the first index
                    if len(refseq_indices) > 1:
                        refseq_entries = [gene_dict[sp][gene][i] for i in refseq_indices]

                        # sort by accession first so that order is not determined by the download order
                        # accessions are sorted numerically based on the digits contained
                        refseq_entries.sort(key=lambda i:
                                            (0, int("".join(findall(r'\d+', i[ORDER_DICT["accession"]]))))
                                            if len(findall(r'\d+', i[ORDER_DICT["accession"]])) > 0
                                            else (1, i[ORDER_DICT["accession"]]))

                        # Then sort by length so that the longest is chosen
                        gene_dict[sp][gene] = sorted(refseq_entries,
                                                     key=lambda i: int(i[ORDER_DICT["length"]]), reverse=True)[0]

                    else:
                        refseq_idx = refseq_indices[0]
                        gene_dict[sp][gene] = gene_dict[sp][gene][refseq_idx]

                # if there aren't any refseq entries, the longest should be the representative sequence
                else:

                    lengths = [int(i[ORDER_DICT["length"]]) for i in gene_dict[sp][gene]]
                    longest_gene = max(lengths)

                    longest_occ = lengths.count(longest_gene)
                    # if all genes are the same length, find the most common sequence from the available sequences
                    if longest_occ > 1:
                        gene_dict[sp][gene] = find_common_seq([i for i in gene_dict[sp][gene] if
                                                               int(i[ORDER_DICT["length"]]) == longest_gene],
                                                              SEQUENCE_DICTIONARY)
                        same_ct += 1
                    # if there's only one longest gene, its the representative sequence
                    else:
                        gene_dict[sp][gene] = sorted(gene_dict[sp][gene],
                                                     key=lambda i: int(i[ORDER_DICT["length"]]), reverse=True)[0]
                        diff_ct += 1
            elif len(gene_dict[sp][gene]) == 1:
                gene_dict[sp][gene] = gene_dict[sp][gene][0]


def find_common_seq(gene_list, seq_dict):
    chosen_gene = []
    record_list = []

    gene_dict = {}
    for gene in gene_list:
        gene_dict["%s.%s" % (gene[ORDER_DICT["accession"]],  gene[ORDER_DICT["gene"]])] = gene
        temp_seq = seq_dict["%s.%s" % (gene[ORDER_DICT["accession"]], gene[ORDER_DICT["gene"]])]
        record_list.append(temp_seq)
    seq_list = [i.seq.count("N") for i in record_list]

    # first, get sequences with the lowest numbers of ambiguous bases
    least_ambig = []
    for i in range(0, len(record_list)):
        if seq_list[i] == min(seq_list):
            least_ambig.append(record_list[i])
    seq_list = [i.seq for i in least_ambig]
    #print(least_ambig)

    # then try to find the most common sequence
    seq_counter = Counter(seq_list).most_common()
    common_count = seq_counter[0][1]
    seq_counter = [i for i in seq_counter if i[1] is common_count]

    # If there is more than one sequence that is the most common, try to get the most recent sequence
    if len(seq_counter) != 1:
        focal_genes = [gene_dict[i.id] for i in least_ambig]
        # Sort by dates so that more recent sequences come first
        old_focal = focal_genes.copy()
        focal_genes.sort(key=lambda i: datetime.strptime(i[ORDER_DICT["date"]], "%d-%b-%Y"), reverse=True)

        if old_focal != focal_genes:
            focal_genes = [i for i in focal_genes if i[ORDER_DICT["date"]] == focal_genes[0][ORDER_DICT["date"]]]

        focal_genes.sort(key=lambda i: (0, int("".join(findall(r'\d+', i[ORDER_DICT["accession"]]))))
                         if len(findall(r'\d+', i[ORDER_DICT["accession"]])) > 0
                         else (1, i[ORDER_DICT["accession"]]))

        # if more than one sequence is the most recent, choose the one at the beginning of the list
        # TODO: this part doesn't have logic behind it except that one sequence has to be chosen!
        chosen_gene = focal_genes[0]
    # if there is a most common sequence, choose the first occurance of that sequence in the list
    else:
        chosen_sequence = seq_counter[0][0]
        for i in least_ambig:
            if i.seq == chosen_sequence:
                chosen_gene = gene_dict[i.id]
                break

    if chosen_gene is []:
        print("FATAL: did not asssign a gene, something went wrong...")
        exit()

    return chosen_gene


def make_acc_matrix(rep_gene_matrix):
    final_acc_matrix = []
    spp_names = list(rep_gene_matrix.keys())
    gene_names = list(rep_gene_matrix[spp_names[0]].keys())

    for i in spp_names:
        temp_list = [i]
        for j in gene_names:
            if not rep_gene_matrix[i][j]:
                temp_list.append("")
            else:
                temp_list.append(rep_gene_matrix[i][j][ORDER_DICT["accession"]])
        final_acc_matrix.append(temp_list)

    gene_names.insert(0, "species")
    final_acc_matrix.insert(0, gene_names)

    return final_acc_matrix


def write_focal_genes(rep_gene_matrix):
    spp_names = list(rep_gene_matrix.keys())
    gene_names = list(rep_gene_matrix[spp_names[0]].keys())
    focal_genes = []
    focal_seqs = []

    for sp in spp_names:
        for gene in gene_names:
            if rep_gene_matrix[sp][gene] != []:
                focal_genes.append("%s.%s" % (rep_gene_matrix[sp][gene][ORDER_DICT["accession"]],
                                              rep_gene_matrix[sp][gene][ORDER_DICT["gene"]]))
    for acc in focal_genes:
        focal_seqs.append(SEQUENCE_DICTIONARY[acc])

    SeqIO.write(focal_seqs, FOCAL_SEQFILE_NAME, "fasta")


def write_2d_matrix(focal_matrix, file_handle):
    with open(file_handle, 'w') as f:
        for i in focal_matrix:
            for j in i:
                f.write(j)
                f.write("\t")
            f.write("\n")


def recalculate_gene_length(gene_list):
    for gene in gene_list:
        gene_seq = SEQUENCE_DICTIONARY["%s.%s" % (gene[ORDER_DICT["accession"]], gene[ORDER_DICT["gene"]])].seq
        new_length = len(gene_seq) - gene_seq.count("N")
        gene[ORDER_DICT["length"]] = new_length


if __name__ == '__main__':
    g_path = Path(args.scraped_genes).resolve()
    g_list = parse_file_nohead(g_path)
    check_gaps(g_list)
    recalculate_gene_length(g_list)

    uniq_genes = get_uniq_element(g_list, "gene")
    uniq_spp = get_uniq_element(g_list, "species")

    aberrant_spp = clean_species(uniq_spp)
    novels, subs = get_weird_spp(uniq_spp)

    uniq_spp = [x for x in uniq_spp if x not in novels and x not in subs]
    # filter out undefined species
    uniq_spp = [x for x in uniq_spp if " sp." not in x]
    # write a taxonomy file
    with open(SPP_FILE_NAME, 'w') as species_file:
        for x in uniq_spp:
            species_file.write(x)
            species_file.write("\n")

    uniq_spp_dict = make_gene_dict(uniq_spp, uniq_genes, g_list)
    # test_sp = "Thorius munificus"
    # print(uniq_spp_dict[test_sp])
    # print("--------------")
    filter_gene_dict(uniq_spp_dict)

    # print(uniq_spp_dict[test_sp])
    write_focal_genes(uniq_spp_dict)
    write_2d_matrix(make_acc_matrix(uniq_spp_dict), "scraped_matrix")

    exit()
    nov_spp_dict = make_gene_dict(novels, uniq_genes, g_list)
    sub_spp_dict = make_gene_dict(subs, uniq_genes, g_list)
