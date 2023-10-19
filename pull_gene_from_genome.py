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

parser = argparse.ArgumentParser(description="Program: Pull Gene from Genome\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("rep_gene_fasta", help='File containing genes from rep_seqs.py', type=str)
parser.add_argument("gene_search_file", help="same file used in gene_scraper.py")
parser.add_argument("blast_db", help='path to a blastn database for the genome of interest', type=str)
parser.add_argument("organism_gtf", help='path to a genome track file for the genome of interest', type=str)
parser.add_argument("--blaste", type=str, default="1e-20",
                    help='path to a blastn database for the genome of interest')

args = parser.parse_args()

rep_genes_seqdict = SeqIO.to_dict(SeqIO.parse(args.rep_gene_fasta, "fasta"))
wiggle = 5000


def run_blastn(db, query):
    print("\nRunning BLASTn with %s as db and %s as query...." % (Path(db).resolve(), Path(query).resolve()))
    # blast run contains only the necessary information
    blastcommand = ['blastn', '-db', str(Path(db).resolve()), "-query", query]

    blastcommand.extend(["-num_threads", "6"])

    blastcommand.extend(["-evalue", args.blaste, "-outfmt", "6"])
    blast = subprocess.run(blastcommand, stdout=subprocess.PIPE, universal_newlines=True)
    blastcommand[len(blastcommand) - 1] = "\"" + blastcommand[len(blastcommand) - 1] + "\""

    print("BLAST was run with these options: ")
    print(" ".join(blastcommand))

    blast = parse_blastn(blast)

    return blast


def parse_blastn(blast_out):
    blast_output = blast_out.stdout.split("\n")
    final_output = []
    for hit in blast_output:
        hit = hit.split("\t")
        for ele in range(0, len(hit)):
            hit[ele] = dv.BLASTN_TYPES[dv.BLASTN_DEFAULT_OUTFMT_SIX[ele]](hit[ele])
        if hit != [""]:
            final_output.append(hit)
            # print(hit)
    return final_output


def parse_gtf(gtf_path):
    gtf_list = utils.parse_file_nohead_tolist(gtf_path)
    for line in gtf_list:
        line[dv.GTF_FORMAT["start"]] = int(line[dv.GTF_FORMAT["start"]])
        line[dv.GTF_FORMAT["end"]] = int(line[dv.GTF_FORMAT["end"]])
        att_dict = {}
        temp_att = line[dv.GTF_FORMAT["attribute"]].split(";")
        for entry in temp_att:
            temp_entry = entry.strip().split(" ", 1)
            if temp_entry != [""]:
                try:
                    att_dict[temp_entry[0]] = temp_entry[1].replace("\"", "")
                except IndexError:
                    print("ERROR: semicolon used inside attribute entry")
                    print(line)
        line[dv.GTF_FORMAT["attribute"]] = att_dict
    return gtf_list


if __name__ == '__main__':
    gtf_table = parse_gtf(args.organism_gtf)
    blast_table = run_blastn(args.blast_db, args.rep_gene_fasta)
    gene_dict = {}
    for k in blast_table:
        try:
            gene_dict[k[0].split(".")[-1]].append(k)
        except KeyError:
            gene_dict.setdefault(k[0].split(".")[-1], [k])

    for gene in gene_dict.keys():
        target_chromosomes = []
        start_stop_values = []
        for hit in gene_dict[gene]:
            # if gene == "12S":
            #    print(hit)
            target_chromosomes.append(hit[1])
        likely_chr = Counter(target_chromosomes).most_common()[0][0]
        for hit in gene_dict[gene]:
            if hit[1] == likely_chr:
                start_stop_values.extend([hit[8], hit[9]])
        focal_gtf_entries = [k for k in gtf_table if k[dv.GTF_FORMAT["seqname"]] == likely_chr and
                             (k[dv.GTF_FORMAT["feature"]].lower() == "gene" or
                              k[dv.GTF_FORMAT["feature"]].lower() == "mrna") and
                             (max(k[dv.GTF_FORMAT["start"]], min(start_stop_values)) <
                              min(k[dv.GTF_FORMAT["end"]], max(start_stop_values))) and
                             (k[dv.GTF_FORMAT["attribute"]]["gene_biotype"].lower() == "protein_coding" or
                              k[dv.GTF_FORMAT["attribute"]]["gene_biotype"].lower() == "rrna")
                             ]
        rrna_entries = [k for k in gtf_table if k[dv.GTF_FORMAT["seqname"]] == likely_chr and
                             k[dv.GTF_FORMAT["feature"]].lower() == "transcript" and
                             (max(k[dv.GTF_FORMAT["start"]], min(start_stop_values)) <
                              min(k[dv.GTF_FORMAT["end"]], max(start_stop_values)))]
        rrna_entries = [k for k in rrna_entries if
                        k[dv.GTF_FORMAT["attribute"]]["transcript_biotype"].lower() == "rrna"]
        print(rrna_entries)
        focal_gtf_entries.extend(rrna_entries)

        found_gene = False
        for line in focal_gtf_entries:
            transcript_q = "gene"
            if line[dv.GTF_FORMAT["feature"]] == "transcript":
                transcript_q = "product"
            desc_in_gois = ""
            gn_in_gois = ""
            try:
                in_gois = dv.MITO_GOIS[line[dv.GTF_FORMAT["attribute"]]["description"].lower()].lower()
            except KeyError:
                pass
            try:
                gn_in_gois = dv.MITO_GOIS[line[dv.GTF_FORMAT["attribute"]][transcript_q].lower()].lower()
            except KeyError:
                pass

            if desc_in_gois == gene.lower() or \
                    gn_in_gois == gene.lower() or \
                    gene.lower() == line[dv.GTF_FORMAT["attribute"]][transcript_q].lower():
                #print("got one!")
                #print(gene)
                #print(line)
                #print("------------")
                gene_dict[gene] = line
                found_gene = True

        if found_gene is False:
            print(gene)
            print(focal_gtf_entries)
            gene_dict[gene] = []

    pprint(gene_dict)