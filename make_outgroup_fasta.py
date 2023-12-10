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
import os
from shutil import rmtree
from pprint import pprint
import shutil

import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Pull Gene from Genome\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("outgroup_genes", type=str,
                    help='path to combined outgroup file made using \'outgroup_pull_genes.py\'')
parser.add_argument("--add_dict", type=str, default="",
                    help="a file containing translations for gene names made using the \'--gbk_dir\' option"
                         " and !!CHECKED FOR INCONSISTENCIES!!")

args = parser.parse_args()

# update GOI dict
if args.add_dict != "":
    utils.update_gene_dict_from_file(args.add_dict, dv.MITO_GOIS)
# make gene description dict
GENE_DESCRIPTION = {gn: [] for gn in sorted(list(set(dv.MITO_GOIS.values())))}
for long_gn in dv.MITO_GOIS.keys():
    GENE_DESCRIPTION[dv.MITO_GOIS[long_gn]].append(long_gn)
for gn in GENE_DESCRIPTION.keys():
    GENE_DESCRIPTION[gn] = sorted(GENE_DESCRIPTION[gn], key=len, reverse=True)[0]

if __name__ == '__main__':
    outgroup_dict = utils.parse_file_head_todict(args.outgroup_genes, 0)
    genomic_parts = {}
    transcript_info = {}
    idx_dict = {k: outgroup_dict["header"][k] for k in range(0, len(outgroup_dict["header"]))}
    for sp in outgroup_dict.keys():
        if sp == "header":
            continue
        for element in range(1, len(outgroup_dict[sp])):
            ele_ref = outgroup_dict[sp][element]
            # If no gene identified, skip
            if ele_ref == "":
                continue
            temp_ele = ele_ref.split(":")
            # if a transcript was identified for the gene, download that instead
            if len(temp_ele) == 1:
                transcript_info[temp_ele[0]] = [idx_dict[element], sp]

            else:
                genomic_parts.setdefault(temp_ele[0], [])

                temp_ele[1] = [int(k) for k in sorted(temp_ele[1].split("-"))]

                if temp_ele[2] == "plus":
                    seq_pol = +1
                else:
                    seq_pol = -1

                element_feature = SeqFeature.FeatureLocation(min(temp_ele[1])-1, max(temp_ele[1]), strand=seq_pol)

                genomic_parts[temp_ele[0]].append([element_feature, idx_dict[element], sp])

    output_dir = "outgroup_chroms"
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    goi_records = []
    print("Fetching transcripts...")
    # fetch transcripts in chunks of 500 sequences at a time
    transcript_chunks = list(utils.chunks(list(transcript_info.keys()), 500))
    chunkfile_names = ["%s/og_transcript_chunk%s.fasta" % (output_dir, k) for k in range(0, len(transcript_chunks))]
    comb_chunkfile = "%s/combined_og_transcript_chunks.fasta" % output_dir

    # only fetch transcripts if the combined chunk file doesn't already exist
    if not Path(comb_chunkfile).exists():
        for acc_chunk in range(0, len(transcript_chunks)):
            handle = Entrez.efetch(db="nucleotide", id=transcript_chunks[acc_chunk], rettype="fasta", retmode="text")
            with open(chunkfile_names[acc_chunk], 'w') as f:
                for k in list(handle):
                    f.write(str(k))

        # combine all of the pulled transcript sequences into one file
        with open(comb_chunkfile, 'wb') as final_file:
            for chunkfile in chunkfile_names:
                with open(chunkfile, 'rb') as temp_file:
                    shutil.copyfileobj(temp_file, final_file)
                os.remove(chunkfile)

    # open and reformat the transcripts downloaded from GenBank
    transcript_seqs = SeqIO.to_dict(SeqIO.parse(comb_chunkfile, "fasta"))
    for tid in transcript_info.keys():
        focal_sp = transcript_info[tid][1]
        focal_gene = transcript_info[tid][0]

        """
        long_desc = "label not included in dictionary"
        try:
            long_desc = GENE_DESCRIPTION[focal_gene]
        except KeyError:
            pass
        """

        new_id = "%s.%s" % (tid, focal_gene)
        new_desc = "%s %s" % (focal_sp, focal_gene)

        transcript_seqs[tid].id = new_id
        transcript_seqs[tid].name = new_id
        transcript_seqs[tid].description = new_desc
        goi_records.append(transcript_seqs[tid])

    loop_iter = 1
    for chrom in genomic_parts.keys():
        if not Path("%s/%s.fasta" % (output_dir, chrom)).exists():
            print("Fetching sequence for %s" % chrom)
            handle = Entrez.efetch(db="nucleotide", id=chrom, rettype="fasta", retmode="text")
            with open("%s/%s.fasta" % (output_dir, chrom), 'w') as f:
                for k in list(handle):
                    f.write(str(k))
        chrom_seq = list(SeqIO.parse("%s/%s.fasta" % (output_dir, chrom), "fasta"))[0].seq

        print("Getting GOIs from %s" % chrom)
        for goi in genomic_parts[chrom]:
            """
            long_desc = "label not included in dictionary"
            try:
                long_desc = GENE_DESCRIPTION[goi[1]]
            except KeyError:
                pass
            """

            goi_records.append(SeqIO.SeqRecord(
                goi[0].extract(chrom_seq),
                id="%s.%s" % (chrom, goi[1]),
                name="%s.%s" % (chrom, goi[1]),
                description="%s %s" % (goi[2], goi[1])
                ))
        print("Done with %s! (%s/%s)\n" % (chrom, loop_iter, len(genomic_parts)))
        loop_iter += 1

    SeqIO.write(goi_records, "combined_outgroups.fasta", "fasta")
