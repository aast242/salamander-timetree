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

parser = argparse.ArgumentParser(description="Program: Everson et al. Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("everson_fasta_dir", help='A *DIRECTORY* containing all of the Everson fastas', type=str)
# this file can be downloaded at the link below:
# https://www.ncbi.nlm.nih.gov/Traces/wgs/?display=download
parser.add_argument("translation_table", type=str,
                    help='A taxonomy translation table from DWW specimen numbers to species')
parser.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")
parser.add_argument("--find_uniq", action='store_true', help="Try to find distinct sequences for each sp")

# >KEBL01000001.1.prdm16 Caecilia tentaculata prdm16

args = parser.parse_args()

if __name__ == '__main__':
    # set up a folder for exports
    fasta_dir_name = Path(args.everson_fasta_dir).stem
    export_dir_name = Path(f"./{fasta_dir_name}_SC-reformat")
    try:
        os.mkdir(export_dir_name)
    except FileExistsError:
        pass

    # make a translation table between DWW IDs and species names
    dww_tt = utils.parse_file_nohead_todict(args.translation_table, 0)
    dww_tt = {i: dww_tt[i][1] for i in dww_tt.keys()}

    # get the names of all of the fasta files to modify
    everson_fastas = [Path(f"{Path(args.everson_fasta_dir)}/{i}") for i in os.listdir(fasta_dir_name)
                      if i.endswith("fna")]

    # iterate through each fasta file
    for file_name in everson_fastas:
        # get the gene name from the fasta file name
        gene_name = file_name.with_suffix("").stem.lower()
        # open the fasta file
        gene_fasta = SeqIO.to_dict(SeqIO.parse(file_name, "fasta"))
        fasta_ids = list(gene_fasta.keys())
        for entry in fasta_ids:
            # turn haplotypes into unique seq IDs
            dww_id = gene_fasta[entry].id.replace("_", ".")
            # try to get the species from the DWW ID using the translation table
            try:
                dww_sp = dww_tt[dww_id.split(".")[0]]
            except KeyError:
                if not args.quiet:
                    print("Key %s not in translation table!" % dww_id.split(".")[0])
                    print("Excluding the following sequence:")
                    print(gene_fasta[entry])
                    print("---------------------------------")
                del gene_fasta[entry]
                continue

            # ungap the sequence
            gene_fasta[entry].seq = gene_fasta[entry].seq.ungap("-")

            # modify the descriptions of the fasta file for SuperCRUNCH
            gene_fasta[entry].id = f"{dww_id}.{gene_name}"
            gene_fasta[entry].name = gene_fasta[entry].id
            gene_fasta[entry].description = f"{dww_sp} {gene_name}"
        gene_fasta = gene_fasta.values()

        # get different sequences for each species
        if args.find_uniq:
            curr_seq_choice = gene_fasta
            p = 0
            while p < 250:
                sp_dict = {}
                seq_dict = {}
                for entry in gene_fasta:
                    entry_sp = " ".join(entry.description.split(" ")[:2])
                    sp_dict.setdefault(entry_sp, [])
                    seq_dict.setdefault(entry_sp, set())

                    sp_dict[entry_sp].append(entry)
                    seq_dict[entry_sp].add(entry.seq)
                uniq_seqs = [list(seq_dict[k]) for k in seq_dict.keys()]
                flat_uniq_seqs = []
                for entry in uniq_seqs:
                    flat_uniq_seqs.extend(entry)
                flat_uniq_seqs = [i[0] for i in list(Counter(flat_uniq_seqs).items()) if i[1] == 1]

                for sp in seq_dict.keys():
                    temp_seqs = seq_dict[sp]
                    keep_seqs = set(i for i in temp_seqs if i in flat_uniq_seqs)
                    if len(keep_seqs) > 0:
                        if len(keep_seqs) > 1:
                            temp_var = list(keep_seqs)[0]
                            keep_seqs = set()
                            keep_seqs.add(temp_var)
                        seq_dict[sp] = keep_seqs

                # find the best set that has the most distinctive sequences
                uniq_seqs = [list(seq_dict[k]) for k in seq_dict.keys()]
                curr_choice = ()
                for i in itertools.product(*uniq_seqs):
                    if len(Counter(i)) > len(curr_choice):
                        curr_choice = i
                curr_choice = list(curr_choice)
                # assign those distinctive sequences to species
                for sp in seq_dict.keys():
                    seq_dict[sp] = list(seq_dict[sp])
                    i = 0
                    while i < len(seq_dict[sp]):
                        # remove the element from the best set that matched and set the species to be only that sequence
                        if seq_dict[sp][i] in curr_choice:
                            curr_choice.remove(seq_dict[sp][i])
                            seq_dict[sp] = [seq_dict[sp][i]]
                            i += 1
                        # never completely empty the list!
                        elif len(seq_dict[sp][i]) == 1:
                            pass
                        # delete things that don't match with anything (shouldn't get here often)
                        else:
                            del seq_dict[sp][i]

                # pull the first matching distinctive sequence from the fasta
                tmp_gene_fasta = []
                for sp in seq_dict.keys():
                    for entry in sp_dict[sp]:
                        if entry.seq == seq_dict[sp][0]:
                            tmp_gene_fasta.append(entry)
                            break

                curr_p_delta_one = [j[1]-1 for j in list(Counter([i.seq for i in tmp_gene_fasta]).items())]
                old_p_delta_one = [j[1]-1 for j in list(Counter([i.seq for i in curr_seq_choice]).items())]
                if sum(curr_p_delta_one) < sum(old_p_delta_one) and \
                        (p > 1 and
                         len([i for i in curr_p_delta_one if i != 0]) >= len([i for i in old_p_delta_one if i != 0])):
                    curr_seq_choice = tmp_gene_fasta
                elif sum(curr_p_delta_one) < sum(old_p_delta_one) and p <= 1:
                    curr_seq_choice = tmp_gene_fasta
                p += 1
            gene_fasta = curr_seq_choice
        #pprint(gene_fasta)

        SeqIO.write(gene_fasta, f"{export_dir_name}/{gene_name}_refmt.fasta", "fasta")

