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
parser.add_argument("blastn_out6", help="", type=str)
parser.add_argument("subj_fasta", help="")
parser.add_argument("query_fasta")
parser.add_argument("sp_name", help="")
parser.add_argument("--bp_bridge", default=200, type=int)

args = parser.parse_args()


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
    blast_output = utils.parse_file_nohead_tolist(blast_out)
    for hit in blast_output:
        for ele in range(0, len(hit)):
            hit[ele] = dv.BLASTN_TYPES[dv.BLASTN_DEFAULT_OUTFMT_SIX[ele]](hit[ele])
    return blast_output


# FROM SUPERCRUNCH, Reference_Blast_Extract.py
# Written by Daniel Portik
# daniel.portik@gmail.com
# January 2019
# Distributed under the
# GNU General Public Lincense
def choose_coord(coords):
    """
    If a coordinate list has multiple non-overlapping intervals,
    find the longest interval and return the single coordinate
    set as a list.
    """
    diff_index = []
    for i, j in enumerate(coords):
        diff = j[-1] - j[0]
        diff_index.append([i, diff])

    diff_index.sort(key=lambda x: x[1], reverse=True)
    selected = coords[(diff_index[0][0])]

    return selected


# FROM SUPERCRUNCH, Reference_Blast_Extract.py
# Written by Daniel Portik
# daniel.portik@gmail.com
# January 2019
# Distributed under the
# GNU General Public Lincense
def merge_coords(interval_list):
    """
    Takes a list of sublists where each sublist is composed of an
    interval [start integer, stop integer]. Assumes that for
    intervals the first number is lower than the second number.
    The sublists are sorted by the start integer, low to
    high, and the search begins with the first interval/sublist.
    If the intervals overlap the starting interval list is
    updated with the new maximum stop integer, and the search
    continues. If there is no overlap, the non-overlapping
    interval becomes the new target for the search, which
    ends after the sublists have been iterated over.

    Example:
    This list of sublists (sorted):
    [ [1,6], [5,15], [7,12], [20,25], [22,35], [35,40] ]

    Creates the merged list:
    [ [1, 6, 15, 15], [20, 25, 35, 40] ]

    Which is cleaned to create the final list:
    [ [1, 15], [20, 40] ]
    """
    # sort the list of lists by their first element (start coordinate)
    interval_list.sort(key=lambda x: x[0])
    # initiate empty list to store updated coordinates
    merged = []
    # iterate over lists
    for interval in interval_list:
        # if merged list is empty, add first interval/sublist
        # otherwise, see if the stop coord of the next interval/sublist
        # is less than the start coord of this (merged) interval
        if not merged or merged[-1][-1] < interval[0]:
            # if there is no overlap, add the non-overlapping interval
            # to the merged list to begin the next round of searching
            merged.append(interval)
        # if there is overlap
        else:
            # find the highest stop coord between the current (merged)
            # interval and the overlapping interval
            int_end = max(merged[-1][-1], interval[-1])
            # append that number to the end of the merged list entry,
            # growing the (merged) interval by the highest stop integer
            # found so far
            merged[-1].append(int_end)
    # take first and last element from each sublist in merged
    # this produces one start and stop coord for every interval
    cleaned = []
    for m in merged:
        cleaned.append([m[0], m[-1]])

    return cleaned


# FROM SUPERCRUNCH, Reference_Blast_Extract.py
# Written by Daniel Portik
# daniel.portik@gmail.com
# January 2019
# Distributed under the
# GNU General Public Lincense
def span_coords(coords, bp_bridge):
    """
    Check if multiple bp intervals are within X bp of each
    other, and if so, merge them and return a new list
    of updated coordinates. X is determined by argument
    bp_bridge, which defaults to 100 unless supplied by user.
    """
    new_coords = []
    for i in range(0, (len(coords) - 1)):
        # check if the end coordinate of the first entry is within X bp
        # of the start coordinate of the subsequent entry
        # if so, merge. X can be set by user (arg --bp_bridge; default is 100)
        if (int(coords[i + 1][0]) - int(coords[i][1])) <= int(bp_bridge):
            updated = [coords[i][0], coords[i + 1][1]]
            new_coords.append(updated)

        else:
            new_coords.append(coords[i])
            new_coords.append(coords[i + 1])

    new_coords.sort(key=lambda x: x[0])

    return new_coords


# FROM SUPERCRUNCH, Reference_Blast_Extract.py
# Written by Daniel Portik
# daniel.portik@gmail.com
# January 2019
# Distributed under the
# GNU General Public Lincense
def one_coord_span(merged_coords, bp_bridge):
    """
    Function for finding the most sensible blast coordinates to use.
    If multiple non-overlapping coordinates are present,
    will first attempt to connect intervals that are separated
    by less than X bp (default 100). If multiple intervals
    persist, then a long stretch of N's is present in the
    sequence (making it low-quality) or it resulted from
    distinct hits from one sequence. In this case, the longest
    interval is retained (which should represent the target locus).
    """
    # NOTE: merged_coords always returns a list of sublists,
    # even if only one sublist, ie [[start,stop]]
    # must return same data structure here
    if len(merged_coords) >= int(2):
        spanned = span_coords(merged_coords, bp_bridge)

        if len(spanned) >= int(2):
            spanned = merge_coords(spanned)

        if len(spanned) >= int(2):
            final_coord = choose_coord(spanned)

        else:
            final_coord = spanned[0]

    else:
        final_coord = merged_coords[0]

    return final_coord


if __name__ == '__main__':
    print("Parsing %s" % args.blastn_out6)
    blast_table = parse_blastn(args.blastn_out6)
    #for k in blast_table[:10]:
    #    print(k)

    query_dict = SeqIO.to_dict(SeqIO.parse(args.query_fasta, "fasta"))
    blast_dict = {}
    for k in blast_table:
        blast_dict.setdefault(k[0].split(".")[-1], [])
        blast_dict[k[0].split(".")[-1]].append(k)

    print("getting coordinate ranges for GOIs")
    gene_list = []
    for k in blast_dict.keys():
        # get the 500 longest hits

        top_one_hundred = sorted(blast_dict[k], key=lambda x: x[3], reverse=True)
        longest_hit = top_one_hundred[0][3]

        top_one_hundred = [x for x in top_one_hundred if x[3] >= (longest_hit * 0.50)]
        half_length_cutoff = int(len(top_one_hundred) * 0.51)
        if half_length_cutoff == 0:
            half_length_cutoff = 1
        # get the 100 top scoring hits from the longest ones
        top_one_hundred = sorted(top_one_hundred, key=lambda x: -x[11])[:half_length_cutoff]
        highest_score = top_one_hundred[0][11]
        temp_id = Counter([x[1] for x in top_one_hundred]).most_common()[0][0]
        # This caused trouble with ND4L
        #for gene in Counter([x[1] for x in top_one_hundred]).keys():
        #    if k in gene.lower():
        #        temp_id = gene

        # get hits that score at least 50% of the higest scoring hit
        high_scoring_hits = [x for x in blast_dict[k] if x[1] == temp_id and x[11] >= (highest_score * 0.5)]
        temp_id_list = [[min([x[8], x[9]]), max([x[8], x[9]])] for x in high_scoring_hits]
        # if the temp list is empty (i.e., the most common hit was not as good as the highest scoring) continue
        if temp_id_list == []:
            print("Skipping %s" % k)
            continue
        # merge all of the coordinates for all blast hits on that chromosome
        merged_coordinates = merge_coords(temp_id_list)
        # get the one longest coordinate range with a span
        merged_coordinates = one_coord_span(merged_coordinates, args.bp_bridge)
        min_hit, max_hit = [], []

        for entry in high_scoring_hits:
            if entry[8] == merged_coordinates[0] or entry[9] == merged_coordinates[0]:
                min_hit = entry
            if entry[8] == merged_coordinates[1] or entry[9] == merged_coordinates[1]:
                max_hit = entry

        # add bases on to the query range so that it matches with the number of bases in the query sequence
        # in this case, the min hit is reversed, meaning we would need to add bases to the end
        # and subtract from the beginning
        min_hit_len = len(query_dict[min_hit[0]].seq)
        start_adjust = min_hit[6] - 1
        end_adjust = min_hit_len - min_hit[7]

        if min_hit[9] < min_hit[8]:
            new_min_hit_range = [min_hit[9] - end_adjust, min_hit[8] + start_adjust]
        else:
            new_min_hit_range = [min_hit[8] - start_adjust, min_hit[9] + end_adjust]

        max_hit_len = len(query_dict[max_hit[0]].seq)
        start_adjust = max_hit[6] - 1
        end_adjust = max_hit_len - max_hit[7]
        if max_hit[9] < max_hit[8]:
            new_max_hit_range = [max_hit[9] - end_adjust, max_hit[8] + start_adjust]
        else:
            new_max_hit_range = [max_hit[8] - start_adjust, max_hit[9] + end_adjust]

        new_comb_range = new_min_hit_range + new_max_hit_range
        new_comb_range = [min(new_comb_range), max(new_comb_range)]
        merged_coordinates = new_comb_range

        # add some wiggle room in
        # wiggle = int((merged_coordinates[1]-merged_coordinates[0]) * 0.1)
        wiggle = 0

        merged_coordinates = [merged_coordinates[0] - wiggle, merged_coordinates[1] + wiggle]

        if merged_coordinates[0] < 0:
            merged_coordinates[0] = 0
        gene_list.append([temp_id, k,
                          SeqFeature.FeatureLocation(merged_coordinates[0], merged_coordinates[1], strand=+1)])

    blast_dict = {}
    for k in gene_list:
        blast_dict.setdefault(k[0], [])
        blast_dict[k[0]].append(k)

    print("Indexing subject fasta. This might take a while...")
    # subject_idx = SeqIO.index(args.subj_fasta, "fasta")
    subject_idx = SeqIO.to_dict(SeqIO.parse(args.subj_fasta, "fasta"))

    print("Pulling gene sequences from fasta...")
    gene_seqs = []
    for subj_contig in blast_dict.keys():
        temp_contig_seq = subject_idx[subj_contig]
        for gene in blast_dict[subj_contig]:
            temp_seq = gene[2].extract(temp_contig_seq).seq
            temp_id = "%s.%s" % (gene[0], gene[1])
            gene_seqs.append(SeqIO.SeqRecord(
                temp_seq,
                id=temp_id,
                name=temp_id,
                description="%s %s" % (args.sp_name, gene[1])
            ))
    gene_seqs.sort(key=lambda x: x.id.split(".")[-1])
    SeqIO.write(gene_seqs, "%s_outgroup-genes.fasta" % args.sp_name.replace(" ", "_"), "fasta")
