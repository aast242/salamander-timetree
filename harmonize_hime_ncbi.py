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
import itertools
import os
import sys
import shutil
from pprint import pprint
from collections import Counter
import time

from defaults import ProgDefaults as dv
import utils

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Salamander Gene Scraper\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("ncbi_matrix", type=str,
                    help='File containing genes to search for')
parser.add_argument("ncbi_fasta", type=str,
                    help="the txid for the taxon you're pulling genes from")
parser.add_argument("ref_outgroup_matrix", type=str,
                    help="")
parser.add_argument("ref_outgroup_fasta", type=str,
                    help="")
parser.add_argument("hime_fasta", type=str,
                    help="performs a search for full mitochondrial genomes")
parser.add_argument("--add_dict", type=str, default="",
                    help="a file containing translations for gene names made using the \'--gbk_dir\' option"
                         " and !!CHECKED FOR INCONSISTENCIES!!")


# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delimited on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()
if args.add_dict != "":
    utils.update_gene_dict_from_file(args.add_dict, dv.MITO_GOIS)

# harmonize nomenclature across datasources
for k in dv.HIME_TO_SYMBOL.keys():
    try:
        dv.HIME_TO_SYMBOL[k] = dv.MITO_GOIS[dv.HIME_TO_SYMBOL[k]]
    except KeyError:
        print("%s wasn't a key!" % k)
SYMBOL_TO_HIME = {dv.HIME_TO_SYMBOL[k]: k for k in dv.HIME_TO_SYMBOL.keys()}


def prep_ncbi_matrix(ncbi_mat_fp):
    hime_gene_symbols = dv.HIME_TO_SYMBOL.values()

    parsed_matrix = utils.parse_file_head_todict(ncbi_mat_fp, 0)
    matrix_header = parsed_matrix["header"]
    del parsed_matrix["header"]
    matrix_header = {i: matrix_header[i] for i in range(0, len(matrix_header))}

    for sp in parsed_matrix.keys():
        parsed_matrix[sp] = {matrix_header[i]: parsed_matrix[sp][i] for i in range(1, len(parsed_matrix[sp]))}
        for symbol in hime_gene_symbols:
            parsed_matrix[sp].setdefault(symbol, '')
        for symbol in sorted(list(set(dv.MITO_GOIS.values()))):
            parsed_matrix[sp].setdefault(symbol, '')

    return parsed_matrix


if __name__ == '__main__':
    ncbi_mat = prep_ncbi_matrix(args.ncbi_matrix)
    ref_outgroup_mat = prep_ncbi_matrix(args.ref_outgroup_matrix)
    ncbi_seqs = SeqIO.to_dict(SeqIO.parse(args.ncbi_fasta, "fasta"))
    hime_seqs = SeqIO.to_dict(SeqIO.parse(args.hime_fasta, "fasta"))
    ref_outgroup_seqs = SeqIO.to_dict(SeqIO.parse(args.ref_outgroup_fasta, "fasta"))

    hime_spp = set(hime_seqs[k].description.replace(k, "").split("HLL_")[0].strip() for k in hime_seqs.keys())
    hime_spp = sorted(list(hime_spp))
    for k in hime_spp:
        if k in ncbi_mat.keys():
            continue
        ncbi_mat[k] = {n: "" for n in ncbi_mat[list(ncbi_mat.keys())[0]].keys()}

    # resolve entries for hime
    for species in ncbi_mat.keys():
        active_hime_seqs = [k for k in hime_seqs.keys() if
                            species == hime_seqs[k].description.replace(k, "").split("HLL_")[0].strip()]
        if active_hime_seqs != []:
            active_hime_seqs = {dv.HIME_TO_SYMBOL[int(k.split(".")[-1].replace("HLL_", ""))]: k
                                for k in active_hime_seqs}
            #pprint(active_hime_seqs)
        if not active_hime_seqs:
            continue

        for gene in active_hime_seqs.keys():
            if ncbi_mat[species][gene] == '':
                try:
                    ncbi_mat[species][gene] = active_hime_seqs[gene]
                except KeyError:
                    pass
            else:
                temp_ncbi = ncbi_seqs[ncbi_mat[species][gene]]
                temp_hime = hime_seqs[active_hime_seqs[gene]]
                # if the sequences are the same length or hime is longer, take the hime seq
                if len(temp_hime.seq) >= len(temp_ncbi.seq):
                    ncbi_mat[species][gene] = active_hime_seqs[gene]

        if active_hime_seqs != {}:
            pass
            #print(species)

    hime_seq_ids = hime_seqs.keys()
    ncbi_seq_ids = ncbi_seqs.keys()

    for k in ref_outgroup_mat.keys():
        if k not in ncbi_mat.keys():
            ncbi_mat[k] = ref_outgroup_mat[k]
            for gene in ncbi_mat[k]:
                if ncbi_mat[k][gene] != "":
                    ncbi_mat[k][gene] = "%s.%s" % (ncbi_mat[k][gene].split(":", 1)[0], gene)
            continue
        for gene in ref_outgroup_mat[k].keys():
            if ".HLL_" in ncbi_mat[k][gene]:
                use_seqs = hime_seqs
            else:
                use_seqs = ncbi_seqs

            acc_num = ref_outgroup_mat[k][gene].split(":", 1)[0]

            if ncbi_mat[k][gene] == "" and ref_outgroup_mat[k][gene] != "":
                ncbi_mat[k][gene] = "%s.%s" % (acc_num, gene)
            elif ncbi_mat[k][gene] == "" and ref_outgroup_mat[k][gene] == "":
                pass
            elif ncbi_mat[k][gene] != "" and ref_outgroup_mat[k][gene] == "":
                pass
            else:
                temp_refout = ref_outgroup_seqs["%s.%s" % (acc_num, gene)]
                temp_ncbi = use_seqs[ncbi_mat[k][gene]]

                if temp_refout.seq == temp_ncbi.seq:
                    pass
                elif len(temp_refout.seq) >= len(temp_ncbi.seq):
                    ncbi_mat[k][gene] = "%s.%s" % (acc_num, gene)

    gene_order = list(ncbi_mat[list(ncbi_mat.keys())[0]].keys())
    gene_order.sort()

    # rename all of the hime markers to have symbolic names
    old_hime_keys = list(hime_seqs.keys())
    for record in old_hime_keys:
        hime_num = int(record.split(".")[-1].replace("HLL_", ""))
        hime_symbol = dv.HIME_TO_SYMBOL[hime_num]

        hime_seqs[record].id = hime_seqs[record].id.replace("HLL_%s" % hime_num, hime_symbol)
        hime_seqs[record].name = hime_seqs[record].name.replace("HLL_%s" % hime_num, hime_symbol)
        hime_seqs[record].description = hime_seqs[record].description.replace("HLL_%s" % hime_num, hime_symbol)
        hime_seqs[hime_seqs[record].id.replace("HLL_%s" % hime_num, hime_symbol)] = hime_seqs[record]
        del hime_seqs[record]

    # combine all of the seq dicts
    ncbi_seqs.update(hime_seqs)
    del hime_seqs
    ncbi_seqs.update(ref_outgroup_seqs)
    del ref_outgroup_seqs

    for k in list(ncbi_seqs.keys()):
        sp = ncbi_seqs[k].description.replace(k, "").replace(k.split(".")[-1], "").strip().replace(" ", "_")
        sp = "%s.%s" % (sp, ncbi_seqs[k].id)
        ncbi_seqs[k].description = ncbi_seqs[k].description.replace(ncbi_seqs[k].id, sp)
        ncbi_seqs[k].id = sp
        ncbi_seqs[k].name = sp

    final_sequences = []
    final_matrix = [["species"]]
    final_matrix[0].extend(gene_order)
    for sp in ncbi_mat.keys():
        temp_entry = [sp]
        for gene in gene_order:
            if ".HLL_" in ncbi_mat[sp][gene]:
                hime_num = int(ncbi_mat[sp][gene].split(".")[-1].replace("HLL_", ""))
                hime_symbol = dv.HIME_TO_SYMBOL[hime_num]
                ncbi_mat[sp][gene] = ncbi_mat[sp][gene].replace("HLL_%s" % hime_num, hime_symbol)
            if ncbi_mat[sp][gene] != "":
                try:
                    final_sequences.append(ncbi_seqs[ncbi_mat[sp][gene]])
                except KeyError:
                    print(sp)
                    print(gene)
                    print(ncbi_mat[sp][gene])
                    exit()
                ncbi_mat[sp][gene] = "%s.%s" % (sp.replace(" ", "_"), ncbi_mat[sp][gene])
            temp_entry.append(ncbi_mat[sp][gene])
        final_matrix.append(temp_entry)

    with open("saltree_final_harmonized_sequences_matrix.txt", "w") as f:
        for line in final_matrix:
            for entry in line:
                f.write("%s\t" % entry)
            f.write("\n")

    SeqIO.write(final_sequences, "saltree_final_harmonized_sequences.fasta", "fasta")
