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
parser.add_argument("final_seqs_preconcat", help='A *DIRECTORY* containing pre-concatenation fastas', type=str)
parser.add_argument("trimmed_seqs_spp_id", type=str, help="A *DIRECTORY* containing the renamed fastas with spp and ID")
parser.add_argument("--starting_seq_dir", type=str, help="A *DIRECTORY* containing the starting sequences for SC",
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


def refmt_aln_file(aln_path, trans_table):
    gene_name = aln_path.with_suffix("").stem.lower().replace("-", "_")
    gene_name = trans_table[gene_name]
    aln_list = utils.parse_file_nohead_tolist(aln_path)
    for line in aln_list:
        if line[0].startswith(">"):
            # get the name of the species from fasta header lines
            sp_name = line[0].split("@")[0].split("_")[0].replace(">", "").strip()
            if sp_name in species_to_change.keys():
                sp_name = species_to_change[sp_name]
            # get the UID of the sequence from the fasta header lines and replace . with _
            sq_uid = line[0].split("@")[1].replace(".", "_").lower()

            # print(f"gene: {gene_name}\nsp: {sp_name}\nseq_uid: {sq_uid}\n\n\n")
            new_seqline = f">{sq_uid}.{gene_name} {sp_name} {gene_name}"
            line[0] = new_seqline
        else:
            # remove spaces and asterisks from alignments
            stripped_seq = line[0].replace(" ", "").replace("*", "")
            line[0] = stripped_seq
    return aln_list


if __name__ == '__main__':
    starting_seq_tt = {}
    if args.starting_seq_dir != "":
        for i in ["%s/%s" % (Path(args.starting_seq_dir), i) for i in list(os.listdir(args.starting_seq_dir))]:

            starting_seq_tt[filename_tt[Path(i).stem]] = sorted(list({record.id: record for record
                                                                      in SeqIO.parse(i, "fasta")}.keys()))
    # set up the translation table
    final_files = ["%s/%s" % (Path(args.final_seqs_preconcat), i)
                   for i in list(os.listdir(args.final_seqs_preconcat))]
    rename_files = ["%s/%s" % (Path(args.trimmed_seqs_spp_id), i)
                    for i in list(os.listdir(args.trimmed_seqs_spp_id))]

    final_files_dict = {}
    for i in final_files:
        spp = sorted(list(SeqIO.to_dict(SeqIO.parse(i, "fasta")).keys()))
        gn = Path(i).stem.split("_oneseq")[0]
        final_files_dict[gn] = spp

    rename_dict = {}
    for i in rename_files:
        gn = Path(i).stem.split("_oneseq")[0]
        spp = sorted(list(SeqIO.to_dict(SeqIO.parse(i, "fasta")).keys()))
        rename_dict[gn] = {}
        for entry in spp:
            tmp = entry.split("_")
            sn = "_".join(tmp[:2])
            uid = "_".join(tmp[2:])
            rename_dict[gn][sn] = uid

    final_info_list = []
    for gene in final_files_dict.keys():
        for sp in final_files_dict[gene]:
            curr_dataset = ""
            try:
                curr_uid = rename_dict[gene][sp]
                curr_uid_minus_gene = rename_dict[gene][sp].replace(".%s" % gene, "")
                for starter_file in starting_seq_tt.keys():
                    if curr_uid in starting_seq_tt[starter_file]:
                        curr_dataset = starter_file
                        break
                final_info_list.append([gene, sp, curr_uid, curr_uid_minus_gene, curr_dataset])
            except KeyError:
                pass

    #pprint(final_info_list[:100])
    with open("seq_info_refmt.tsv", "w") as f:
        for line in final_info_list:
            for ele in line:
                f.write(ele)
                f.write("\t")
            f.write("\n")