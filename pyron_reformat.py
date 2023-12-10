#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez, SeqIO, SeqRecord, Seq, SeqFeature
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

ACCEPTABLE_SUBPARSERS = ["2020", "2022", "rename"]
# ------------------------------------------------------------------------ #
# initialize the top-level parser and add the subparsers to the top parser #
# ------------------------------------------------------------------------ #
top_parser = argparse.ArgumentParser('pyron_parser', add_help=False)

# print help message if prog is called w/o subparser OR if subparser doesn't match any of the defined ones #
if len(top_parser.parse_known_args()[1]) == 0 or \
        top_parser.parse_known_args()[1][0] not in ACCEPTABLE_SUBPARSERS:
    print("Usage: pyron_reformat.py 2020|2022|rename")
    exit()

subparsers = top_parser.add_subparsers(dest='subparser_id')

# ----------------------------------------- #
# Define the subparser for the 2020 dataset #
# ----------------------------------------- #
parser_2020 = subparsers.add_parser('2020',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description='Reformat the Pyron et al. 2020 dataset',
                                    usage='pyron_reformat.py 2020 <ahe_phylip> <ahe_parts> '
                                          '<mito_phylip> <mito_parts> <translation_table> [options]')

parser_2020.add_argument("ahe_phylip", type=str,
                         help='A phylip file containing the AHE loci (desmo161_381loci_trim_names.phy)')
parser_2020.add_argument("ahe_parts", type=str,
                         help='A file containing the partitions for the AHE loci '
                              '(desmo161_381loci_trim_names.charSets)')
parser_2020.add_argument("mito_phylip", type=str,
                         help='A phylip file containing the mitochondrial loci (desmo161_mt.phy)')
parser_2020.add_argument("mito_parts", type=str,
                         help='A file containing the partitions for the mitochondrial loci (desmo161_mt.charSets)')
parser_2020.add_argument("translation_table", type=str,
                         help='A taxonomy translation table to latin names')
parser_2020.add_argument("--exclude", type=str, default="",
                         help="file with gene_ids to exclude")
parser_2020.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")

# ----------------------------------------- #
# Define the subparser for the 2022 dataset #
# ----------------------------------------- #
parser_2022 = subparsers.add_parser('2022',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description='Reformat the Pyron et al. 2022 dataset',
                                    usage='pyron_reformat.py 2022 <ahe_phylip_dir> <translation_table> [options]')

parser_2022.add_argument("ahe_phylip_dir", type=str,
                         help='Path to directory containing individual phylip files for each AHE locus')
parser_2022.add_argument("translation_table", type=str, help='A taxonomy translation table to latin names')
parser_2022.add_argument("--quiet", action="store_true", help="Don't print when keys aren't present")
parser_2022.add_argument("--exclude", type=str, default="",
                         help="file with gene_ids to exclude")

parser_rename = subparsers.add_parser('rename',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description='Reformat the Pyron et al. 2022 dataset',
                                    usage='pyron_reformat.py 2022 <ahe_phylip_dir> <translation_table> [options]')

parser_rename.add_argument("comb_fasta", type=str)
parser_rename.add_argument("species_rename", type=str)
parser_rename.add_argument("ahe_rename", type=str)
parser_rename.add_argument("ahe_to_symbol", type=str)

# ------------------------------ #
# parse the inputs from the user #
# ------------------------------ #
args = top_parser.parse_args()

# >KEBL01000001.1.prdm16 Caecilia tentaculata prdm16

# translation table for homologous genes to Everson et al. 2021 determined via blast
PYRON_TRANSTABLE = {
    'desmo900-l1': 'desmo2020-l234',
    'desmo900-l10': 'desmo2020-l161',
    'desmo900-l100': 'desmo2020-l96',
    'desmo900-l101': 'desmo2020-l18',
    'desmo900-l102': 'desmo2020-l272',
    'desmo900-l103': 'desmo2020-l317',
    'desmo900-l104': 'desmo2020-l150',
    'desmo900-l105': 'desmo2020-l139',
    'desmo900-l106': 'desmo2020-l141',
    'desmo900-l107': 'desmo2020-l340',
    'desmo900-l108': 'desmo2020-l215',
    'desmo900-l109': 'desmo2020-l299',
    'desmo900-l11': 'desmo2020-l260',
    'desmo900-l110': 'desmo2020-l118',
    'desmo900-l111': 'desmo2020-l47',
    'desmo900-l112': 'desmo2020-l256',
    'desmo900-l113': 'desmo2020-l5',
    'desmo900-l114': 'desmo2020-l123',
    'desmo900-l115': 'desmo2020-l239',
    'desmo900-l116': 'desmo2020-l116',
    'desmo900-l117': 'desmo2020-l199',
    'desmo900-l118': 'desmo2020-l223',
    'desmo900-l119': 'desmo2020-l304',
    'desmo900-l12': 'desmo2020-l229',
    'desmo900-l120': 'desmo2020-l175',
    'desmo900-l121': 'desmo2020-l273',
    'desmo900-l122': 'desmo2020-l271',
    'desmo900-l123': 'desmo2020-l267',
    'desmo900-l124': 'desmo2020-l1',
    'desmo900-l125': 'desmo2020-l298',
    'desmo900-l126': 'desmo2020-l301',
    'desmo900-l127': 'desmo2020-l98',
    'desmo900-l128': 'desmo2020-l111',
    'desmo900-l129': 'desmo2020-l144',
    'desmo900-l13': 'desmo2020-l220',
    'desmo900-l130': 'desmo2020-l288',
    'desmo900-l131': 'desmo2020-l124',
    'desmo900-l132': 'desmo2020-l152',
    'desmo900-l133': 'desmo2020-l329',
    'desmo900-l134': 'desmo2020-l154',
    'desmo900-l135': 'desmo2020-l224',
    'desmo900-l136': 'desmo2020-l186',
    'desmo900-l137': 'desmo2020-l126',
    'desmo900-l138': 'desmo2020-l41',
    'desmo900-l139': 'desmo2020-l351',
    'desmo900-l14': 'desmo2020-l107',
    'desmo900-l140': 'desmo2020-l274',
    'desmo900-l141': 'desmo2020-l357',
    'desmo900-l142': 'desmo2020-l172',
    'desmo900-l143': 'desmo2020-l202',
    'desmo900-l145': 'desmo2020-l127',
    'desmo900-l146': 'desmo2020-l151',
    'desmo900-l147': 'desmo2020-l62',
    'desmo900-l148': 'desmo2020-l131',
    'desmo900-l149': 'desmo2020-l338',
    'desmo900-l15': 'desmo2020-l29',
    'desmo900-l150': 'desmo2020-l176',
    'desmo900-l151': 'desmo2020-l168',
    'desmo900-l152': 'desmo2020-l349',
    'desmo900-l153': 'desmo2020-l27',
    'desmo900-l154': 'desmo2020-l82',
    'desmo900-l155': 'desmo2020-l307',
    'desmo900-l156': 'desmo2020-l295',
    'desmo900-l157': 'desmo2020-l17',
    'desmo900-l158': 'desmo2020-l56',
    'desmo900-l159': 'desmo2020-l330',
    'desmo900-l160': 'desmo2020-l180',
    'desmo900-l161': 'desmo2020-l67',
    'desmo900-l162': 'desmo2020-l110',
    'desmo900-l163': 'desmo2020-l289',
    'desmo900-l164': 'desmo2020-l264',
    'desmo900-l165': 'desmo2020-l97',
    'desmo900-l166': 'desmo2020-l328',
    'desmo900-l167': 'desmo2020-l249',
    'desmo900-l168': 'desmo2020-l177',
    'desmo900-l169': 'desmo2020-l233',
    'desmo900-l17': 'desmo2020-l266',
    'desmo900-l170': 'desmo2020-l59',
    'desmo900-l171': 'desmo2020-l230',
    'desmo900-l172': 'desmo2020-l80',
    'desmo900-l173': 'desmo2020-l90',
    'desmo900-l174': 'desmo2020-l138',
    'desmo900-l175': 'desmo2020-l106',
    'desmo900-l176': 'desmo2020-l137',
    'desmo900-l177': 'desmo2020-l269',
    'desmo900-l178': 'desmo2020-l79',
    'desmo900-l179': 'desmo2020-l125',
    'desmo900-l18': 'desmo2020-l166',
    'desmo900-l180': 'desmo2020-l66',
    'desmo900-l182': 'desmo2020-l237',
    'desmo900-l183': 'desmo2020-l109',
    'desmo900-l184': 'desmo2020-l33',
    'desmo900-l185': 'desmo2020-l19',
    'desmo900-l186': 'desmo2020-l204',
    'desmo900-l187': 'desmo2020-l25',
    'desmo900-l188': 'desmo2020-l198',
    'desmo900-l189': 'desmo2020-l184',
    'desmo900-l19': 'desmo2020-l76',
    'desmo900-l190': 'desmo2020-l30',
    'desmo900-l191': 'desmo2020-l294',
    'desmo900-l192': 'desmo2020-l214',
    'desmo900-l193': 'desmo2020-l34',
    'desmo900-l194': 'desmo2020-l346',
    'desmo900-l195': 'desmo2020-l57',
    'desmo900-l196': 'desmo2020-l276',
    'desmo900-l197': 'desmo2020-l38',
    'desmo900-l198': 'desmo2020-l227',
    'desmo900-l199': 'desmo2020-l49',
    'desmo900-l2': 'desmo2020-l136',
    'desmo900-l20': 'desmo2020-l319',
    'desmo900-l200': 'desmo2020-l20',
    'desmo900-l201': 'desmo2020-l100',
    'desmo900-l203': 'desmo2020-l362',
    'desmo900-l204': 'desmo2020-l210',
    'desmo900-l205': 'desmo2020-l149',
    'desmo900-l206': 'desmo2020-l284',
    'desmo900-l207': 'desmo2020-l121',
    'desmo900-l208': 'desmo2020-l15',
    'desmo900-l209': 'desmo2020-l119',
    'desmo900-l21': 'desmo2020-l65',
    'desmo900-l210': 'desmo2020-l265',
    'desmo900-l211': 'desmo2020-l142',
    'desmo900-l212': 'desmo2020-l88',
    'desmo900-l213': 'desmo2020-l93',
    'desmo900-l214': 'desmo2020-l45',
    'desmo900-l215': 'desmo2020-l146',
    'desmo900-l216': 'desmo2020-l165',
    'desmo900-l217': 'desmo2020-l259',
    'desmo900-l218': 'desmo2020-l312',
    'desmo900-l219': 'desmo2020-l61',
    'desmo900-l22': 'desmo2020-l240',
    'desmo900-l220': 'desmo2020-l318',
    'desmo900-l221': 'desmo2020-l354',
    'desmo900-l222': 'desmo2020-l339',
    'desmo900-l223': 'desmo2020-l187',
    'desmo900-l224': 'desmo2020-l310',
    'desmo900-l225': 'desmo2020-l179',
    'desmo900-l226': 'desmo2020-l197',
    'desmo900-l227': 'desmo2020-l40',
    'desmo900-l228': 'desmo2020-l306',
    'desmo900-l23': 'desmo2020-l296',
    'desmo900-l230': 'desmo2020-l163',
    'desmo900-l231': 'desmo2020-l37',
    'desmo900-l232': 'desmo2020-l9',
    'desmo900-l233': 'desmo2020-l78',
    'desmo900-l234': 'desmo2020-l54',
    'desmo900-l235': 'desmo2020-l212',
    'desmo900-l24': 'desmo2020-l348',
    'desmo900-l25': 'desmo2020-l171',
    'desmo900-l26': 'desmo2020-l278',
    'desmo900-l27': 'desmo2020-l86',
    'desmo900-l28': 'desmo2020-l194',
    'desmo900-l29': 'desmo2020-l46',
    'desmo900-l3': 'desmo2020-l302',
    'desmo900-l30': 'desmo2020-l316',
    'desmo900-l31': 'desmo2020-l344',
    'desmo900-l32': 'desmo2020-l360',
    'desmo900-l33': 'desmo2020-l102',
    'desmo900-l34': 'desmo2020-l281',
    'desmo900-l35': 'desmo2020-l28',
    'desmo900-l36': 'desmo2020-l336',
    'desmo900-l37': 'desmo2020-l192',
    'desmo900-l38': 'desmo2020-l155',
    'desmo900-l39': 'desmo2020-l213',
    'desmo900-l4': 'desmo2020-l71',
    'desmo900-l40': 'desmo2020-l101',
    'desmo900-l41': 'desmo2020-l280',
    'desmo900-l42': 'desmo2020-l326',
    'desmo900-l43': 'desmo2020-l358',
    'desmo900-l44': 'desmo2020-l217',
    'desmo900-l45': 'desmo2020-l170',
    'desmo900-l46': 'desmo2020-l231',
    'desmo900-l47': 'desmo2020-l200',
    'desmo900-l48': 'desmo2020-l173',
    'desmo900-l49': 'desmo2020-l261',
    'desmo900-l5': 'desmo2020-l191',
    'desmo900-l50': 'desmo2020-l282',
    'desmo900-l51': 'desmo2020-l55',
    'desmo900-l52': 'desmo2020-l64',
    'desmo900-l53': 'desmo2020-l2',
    'desmo900-l54': 'desmo2020-l209',
    'desmo900-l55': 'desmo2020-l3',
    'desmo900-l56': 'desmo2020-l350',
    'desmo900-l57': 'desmo2020-l178',
    'desmo900-l58': 'desmo2020-l50',
    'desmo900-l59': 'desmo2020-l232',
    'desmo900-l6': 'desmo2020-l147',
    'desmo900-l60': 'desmo2020-l315',
    'desmo900-l61': 'desmo2020-l253',
    'desmo900-l62': 'desmo2020-l238',
    'desmo900-l63': 'desmo2020-l134',
    'desmo900-l64': 'desmo2020-l270',
    'desmo900-l65': 'desmo2020-l215',
    'desmo900-l66': 'desmo2020-l277',
    'desmo900-l67': 'desmo2020-l52',
    'desmo900-l68': 'desmo2020-l263',
    'desmo900-l69': 'desmo2020-l189',
    'desmo900-l7': 'desmo2020-l255',
    'desmo900-l70': 'desmo2020-l193',
    'desmo900-l71': 'desmo2020-l183',
    'desmo900-l72': 'desmo2020-l331',
    'desmo900-l73': 'desmo2020-l244',
    'desmo900-l74': 'desmo2020-l140',
    'desmo900-l75': 'desmo2020-l60',
    'desmo900-l76': 'desmo2020-l8',
    'desmo900-l77': 'desmo2020-l241',
    'desmo900-l78': 'desmo2020-l188',
    'desmo900-l79': 'desmo2020-l84',
    'desmo900-l8': 'desmo2020-l347',
    'desmo900-l80': 'desmo2020-l108',
    'desmo900-l81': 'desmo2020-l85',
    'desmo900-l82': 'desmo2020-l130',
    'desmo900-l83': 'desmo2020-l314',
    'desmo900-l84': 'desmo2020-l337',
    'desmo900-l85': 'desmo2020-l132',
    'desmo900-l86': 'desmo2020-l246',
    'desmo900-l87': 'desmo2020-l290',
    'desmo900-l88': 'desmo2020-l293',
    'desmo900-l9': 'desmo2020-l133',
    'desmo900-l90': 'desmo2020-l218',
    'desmo900-l91': 'desmo2020-l283',
    'desmo900-l92': 'desmo2020-l228',
    'desmo900-l93': 'desmo2020-l300',
    'desmo900-l94': 'desmo2020-l53',
    'desmo900-l95': 'desmo2020-l275',
    'desmo900-l96': 'desmo2020-l226',
    'desmo900-l97': 'desmo2020-l219',
    'desmo900-l98': 'desmo2020-l195',
    'desmo900-l99': 'desmo2020-l355'
}


def parse_charset(charset_path):
    parse_list = {}
    cycle = 0
    with open(charset_path, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = [k.strip().lower() for k in line.rstrip().split("=")]
            # remove one so that things extract properly from sequences down the line:
            # SeqFeature uses 0 based indexing
            feat_range = [int(k) - 1 for k in line[1].split("-")]
            parse_list[line[0]] = SeqFeature.FeatureLocation(feat_range[0], feat_range[1])
            cycle += 1
    return parse_list


# default SeqIO doesn't like the way the phylips are formatted
def parse_phylip(phylip_path):
    parse_list = []
    cycle = 0
    with open(phylip_path, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue
            if cycle != 0:
                line = line.rstrip().split()
                parse_list.append(SeqRecord.SeqRecord(seq=Seq.Seq(line[1]),
                                                      id=line[0],
                                                      description=""))
            cycle += 1
    return parse_list


def resolve_taxonomy(seq_list, tt, name_part, add_dab):
    sp_seq = 0
    if args.exclude != "":
        exclude_ids = utils.parse_file_nohead_tolist(args.exclude)
    while sp_seq < len(seq_list):
        dab_id = seq_list[sp_seq].id.split("_")
        if not dab_id[name_part] == "P":
            if add_dab:
                dab_id = "DAB%s" % dab_id[name_part]
            else:
                dab_id = dab_id[name_part]
        else:
            dab_id = "_".join(dab_id[name_part:name_part + 2])

        try:
            new_sp_id = tt[dab_id][4].replace("_", "-")
            if new_sp_id == "hubrichti" and dab_id.startswith("P"):
                new_sp_id = "Phaeognathus %s" % new_sp_id
            else:
                new_sp_id = "Desmognathus %s" % new_sp_id

            seq_list[sp_seq].id = dab_id
            seq_list[sp_seq].description = new_sp_id

            sp_seq += 1
        except KeyError:
            if not args.quiet:
                print("----------------")
                print(dab_id)
                print(seq_list[sp_seq])
                print("----------------")
            del seq_list[sp_seq]


def extract_genes(partition_feats, seq_list):
    # extract genes
    gene_extracts = {}
    for name in partition_feats.keys():
        gene_extracts.setdefault(name, [])
        for seq_iso in seq_list:
            curr_extract = partition_feats[name].extract(seq_iso)
            curr_extract.id = "%s.%s" % (curr_extract.id, name)
            curr_extract.name = curr_extract.id
            curr_extract.description = "%s %s" % (curr_extract.description, name)

            gene_extracts[name].append(curr_extract)
    return gene_extracts


def resolve_redundancy(seq_list):
    print("initial list length: %s" % len(seq_list))
    final_seqs = []
    seq_names = [k.id for k in seq_list]
    seq_counts = Counter(seq_names)
    # get the sequences that have multiple entries
    uniq_seqs = [k for k in seq_counts.keys() if seq_counts[k] == 1]
    seq_counts = {k: 0 for k in seq_counts.keys() if seq_counts[k] > 1}
    print("number of redundant sequences: %s" % len(seq_counts.keys()))
    # get the longest sequence for each redundant sequence
    for seq_entry in seq_list:
        try:
            curr_length = len(seq_entry.seq)
            if curr_length > seq_counts[seq_entry.id]:
                seq_counts[seq_entry.id] = curr_length
        except KeyError:
            pass
    print("finished finding the longest sequences")
    redundant_seqids = list(seq_counts.keys())
    count_var = 0
    for seq_entry in seq_list:
        if count_var % 10000 == 0:
            print(count_var)
        # if the sequence is unique, just add it to the final list
        if seq_entry.id in uniq_seqs:
            final_seqs.append(seq_entry)
        # if the sequence has duplicates
        elif seq_entry.id in redundant_seqids:
            # add the longest sequence
            if len(seq_entry.seq) == seq_counts[seq_entry.id]:
                final_seqs.append(seq_entry)
                # remove that sequence id from the list of acceptables
                redundant_seqids.remove(seq_entry.id)
        count_var += 1
    print("finished eliminating redundancies")
    print("final list length: %s" % len(final_seqs))

    return final_seqs


if __name__ == '__main__':
    if args.subparser_id == "rename":
        messy_desmog = list(SeqIO.parse(args.comb_fasta, "fasta"))
        sp_tt = {i[0]: i[1] for i in utils.parse_file_nohead_tolist(args.species_rename)}
        gene_tt = {i[0]: i[1] for i in utils.parse_file_nohead_tolist(args.ahe_rename)}
        symbol_tt = {i[0]: i[1] for i in utils.parse_file_nohead_tolist(args.ahe_to_symbol)}

        sp_keys = sp_tt.keys()
        gene_keys = gene_tt.keys()
        symbol_keys = symbol_tt.keys()

        for entry in messy_desmog:
            space_split_entry = entry.description.split(" ")
            curr_gene = space_split_entry[-1]
            curr_specific = space_split_entry[-2]

            if curr_gene in gene_keys:
                entry.id = entry.id.replace(curr_gene, gene_tt[curr_gene])
                entry.name = entry.name.replace(curr_gene, gene_tt[curr_gene])
                entry.description = entry.description.replace(curr_gene, gene_tt[curr_gene])
                # set curr gene to the renamed one so that translation to symbol works properly
                curr_gene = gene_tt[curr_gene]
            if curr_specific in sp_keys:
                entry.id = entry.id.replace(curr_specific, sp_tt[curr_specific])
                entry.name = entry.name.replace(curr_specific, sp_tt[curr_specific])
                entry.description = entry.description.replace(curr_specific, sp_tt[curr_specific])
            if curr_gene in symbol_keys:
                entry.id = entry.id.replace(curr_gene, symbol_tt[curr_gene])
                entry.name = entry.name.replace(curr_gene, symbol_tt[curr_gene])
                entry.description = entry.description.replace(curr_gene, symbol_tt[curr_gene])
        messy_desmog = resolve_redundancy(messy_desmog)
        SeqIO.write(messy_desmog, "%s.rename" % args.comb_fasta, "fasta")
        exit()

    if args.exclude != "":
        exclude_ids = utils.parse_file_nohead_tolist(args.exclude)
        exclude_ids = [i[0] for i in exclude_ids]

    # make a translation table for species names
    desmog_tt = utils.parse_file_head_todict(args.translation_table, 2)

    # dummy variables so pycharm stops yelling
    analysis_files = {}
    mito_seqs = {}

    if args.subparser_id == "2020":
        # set up a folder for exports
        export_dir_name = Path("./pyron2020_SC-reformat")

        # get mitochondrial genes and partitions
        mito_partitions = parse_charset(args.mito_parts)
        mito_seqs = parse_phylip(args.mito_phylip)
        # update taxonomy
        resolve_taxonomy(mito_seqs, desmog_tt, 0, True)
        # extract genes
        mito_seqs = extract_genes(mito_partitions, mito_seqs)

        nucl_partitions = parse_charset(args.ahe_parts)
        nucl_seqs = parse_phylip(args.ahe_phylip)
        # update taxonomy
        resolve_taxonomy(nucl_seqs, desmog_tt, 1, True)
        # extract genes
        nucl_seqs = extract_genes(nucl_partitions, nucl_seqs)

        # combine the two datasets into one
        mito_seqs.update(nucl_seqs)

    else:  # args.subparser_id == "2022":
        # set up a folder for exports
        export_dir_name = Path("./pyron2022_SC-reformat")
        # parse each of the phylip files
        analysis_files = {Path(i).stem.replace("_", "-").lower(): parse_phylip("%s/%s" % (Path(args.ahe_phylip_dir), i))
                          for i in os.listdir(Path(args.ahe_phylip_dir)) if i.endswith(".phy")}
        for locus_name in analysis_files.keys():
            # change names
            resolve_taxonomy(analysis_files[locus_name], desmog_tt, 1, False)
            # update fasta format
            for entry in analysis_files[locus_name]:
                entry.id = "%s.%s" % (entry.id, locus_name)
                entry.name = entry.id
                entry.description = "%s %s" % (entry.description, locus_name)

    # make the export directory
    try:
        os.mkdir(export_dir_name)
    except FileExistsError:
        pass

    if args.subparser_id == "2020":

        for gene_name in mito_seqs.keys():
            i = 0
            while i < len(mito_seqs[gene_name]):
                mito_seqs[gene_name][i].seq = mito_seqs[gene_name][i].seq.ungap("-")
                if len(mito_seqs[gene_name][i].seq) < 100:
                    del mito_seqs[gene_name][i]
                    i -= 1
                elif args.exclude != "":
                    if mito_seqs[gene_name][i].id.split(".")[-1] in exclude_ids:
                        print("excluding %s" % mito_seqs[gene_name][i].id.split(".")[-1])
                        mito_seqs[gene_name] = []
                i += 1
            if len(mito_seqs[gene_name]) > 0:
                SeqIO.write(mito_seqs[gene_name], f"{export_dir_name}/{gene_name}_refmt.fasta", "fasta")
    else:  # args.subparser_id == "2022":
        for gene_name in analysis_files.keys():
            i = 0
            while i < len(analysis_files[gene_name]):
                analysis_files[gene_name][i].seq = analysis_files[gene_name][i].seq.ungap("-")
                if len(analysis_files[gene_name][i].seq) < 100:
                    del analysis_files[gene_name][i]
                    i -= 1
                i += 1
            SeqIO.write(analysis_files[gene_name], f"{export_dir_name}/{gene_name}_refmt.fasta", "fasta")

# TODO: make exclude work