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

parser = argparse.ArgumentParser(description="Program: Prune Outgroup Seqs\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("fasta_dir", help='A directory containing single-copy seqs for genes', type=str)
args = parser.parse_args()
ambig_genes = ["desmo2020-l", "desmo900-l"]
outgroup_species = ["Latimeria_chalumnae", "Anolis_carolinensis", "Homo_sapiens", "Chrysemys_picta", "Gallus_gallus",
                    "Rhinatrema_bivittatum", "Geotrypetes_seraphini", "Microcaecilia_unicolor",
                    "Xenopus_tropicalis"]

if __name__ == '__main__':
    fasta_list = os.listdir(args.fasta_dir)
    fasta_prefix = str(Path(args.fasta_dir).resolve())
    print(fasta_prefix)
    exclude_ids = []
    for gene_file in fasta_list:
        seq_recs = list(SeqIO.parse("%s/%s" % (fasta_prefix, gene_file), "fasta"))
        outgroup_seqs = []
        for rec in seq_recs:
            curr_sp = "_".join(rec.id.split("_")[:2])
            if curr_sp in outgroup_species:
                outgroup_seqs.append(rec)
        for i in range(0, len(outgroup_seqs)):
            curr_sp = "_".join(outgroup_seqs[i].id.split("_")[:2])
            gene_name = outgroup_seqs[i].id.replace("%s_" % curr_sp, "")
            gene_name = ".".join(gene_name.split(".")[:-1])
            if len(gene_name.split("_")) != 4:
                locus_name = ""
            else:
                locus_name = gene_name.split("_")[2].lower()
            outgroup_seqs[i] = [outgroup_seqs[i].id, curr_sp, gene_name, locus_name]
        current_gene = gene_file.split("_")[0]
        if current_gene == "ahe":
            current_gene = gene_file.split("_")[1]

        outgroup_counter = dict(Counter([i[3] for i in outgroup_seqs if i[3] != ""]))
        try:
            if len(outgroup_counter) > 0:
                outgroup_match_ct = outgroup_counter[current_gene]
                if outgroup_match_ct < 4:
                    print("only %s outgroups matched %s..." % (outgroup_match_ct, current_gene))
                    print(outgroup_counter)
                    potentials = list(outgroup_counter.keys())
                    potentials.remove(current_gene)
                    for potential_id in potentials:
                        syn_q = input("is %s a valid synonym? (y/n): " % potential_id).lower()
                        if syn_q.startswith("y"):
                            pass
                        else:
                            for j in outgroup_seqs:
                                if j[3] == potential_id:
                                    exclude_ids.append(j[0].replace("%s_" % j[1], ""))
                    print()
        except KeyError:
            # If the gene starts with desmo, we don't know what it is.
            # Require all outgroups to match or be synonyms to keep any outgroup sequences
            if current_gene.startswith("desmo"):
                print(current_gene)
                print(outgroup_counter)
                if len(outgroup_counter) > 1:
                    syn_q = input("are all of the outgroup sequences synonymous? (y/n): ").lower()
                    if syn_q.startswith("y"):
                        pass
                    else:
                        for j in outgroup_seqs:
                            exclude_ids.append(j[0].replace("%s_" % j[1], ""))

                print()
            else:
                print("outgroup seqs, but none match %s" % current_gene)
                print(outgroup_counter)
                for potential_id in outgroup_counter.keys():
                    syn_q = input("is %s a valid synonym? (y/n): " % potential_id).lower()
                    if syn_q.startswith("y"):
                        pass
                    else:
                        temp_exclude = [j[0].replace("%s_" % j[1], "") for j in outgroup_seqs if j[3] == potential_id]
                        print(temp_exclude)
                        exclude_ids.extend(temp_exclude)
                print()

    with open("outgroup_pruning_excludes.txt", "w") as f:
        for line in exclude_ids:
            f.write(line)
            f.write("\n")
