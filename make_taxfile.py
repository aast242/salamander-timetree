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
import os
import itertools
import sys
from shutil import rmtree
from pprint import pprint

import utils
from defaults import ProgDefaults as dv

Entrez.email = dv.ENTREZ_EMAIL
Entrez.api_key = dv.ENTREZ_API

parser = argparse.ArgumentParser(description="Program: Salamander Taxfile\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("scraped_genes_fasta", help='File containing genes from gene_scraper.py', type=str)
parser.add_argument("--pre_fmt", action='store_true',
                    help="Assumes input FASTA is preformatted for SuperCRUNCH and doesn't reformat")
parser.add_argument("--collapse_subs", action='store_true',
                    help="Collapses subspecies into binomal names, useful if not analysing subspecies")
parser.add_argument("--include_odds", action='store_true',
                    help="includes undefined species (e.g., \'sp.\', \'aff.\', etc.) in the final output")
parser.add_argument("--exclude", type=str, default="",
                    help="Remove specified loci from the input file")

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delieted on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()


def get_uniq_element(parsed_file, element):
    return sorted(list(set(i[dv.ORDER_DICT[element]] for i in parsed_file)))


def clean_species(spp_list):
    removed_spp = []
    i = 0
    while i < len(spp_list):
        #print([j in spp_list[i] for j in EXCLUDE_TERMS])
        if any([j in spp_list[i] for j in dv.EXCLUDE_TERMS]):
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
        elif args.pre_fmt:
            if len(entry.split(" ")[-1].split("_")) > 1:
                novel_spp.append(entry)
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
            final_dict[entry[dv.ORDER_DICT["species"]]][entry[dv.ORDER_DICT["gene"]]].append(entry)
        except KeyError:
            pass
    return final_dict


def find_common_seq(gene_list, seq_dict):
    chosen_gene = []
    record_list = []

    gene_dict = {}
    for gene in gene_list:
        gene_dict["%s.%s" % (gene[dv.ORDER_DICT["accession"]],  gene[dv.ORDER_DICT["gene"]])] = gene
        temp_seq = seq_dict["%s.%s" % (gene[dv.ORDER_DICT["accession"]], gene[dv.ORDER_DICT["gene"]])]
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
        focal_genes.sort(key=lambda i: datetime.strptime(i[dv.ORDER_DICT["date"]], "%d-%b-%Y"), reverse=True)

        if old_focal != focal_genes:
            focal_genes = [i for i in focal_genes if i[dv.ORDER_DICT["date"]] == focal_genes[0][dv.ORDER_DICT["date"]]]

        focal_genes.sort(key=lambda i: (0, int("".join(findall(r'\d+', i[dv.ORDER_DICT["accession"]]))))
                         if len(findall(r'\d+', i[dv.ORDER_DICT["accession"]])) > 0
                         else (1, i[dv.ORDER_DICT["accession"]]))

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


def write_2d_matrix(focal_matrix, file_handle):
    with open(file_handle, 'w') as f:
        for i in focal_matrix:
            for j in i:
                f.write(j)
                f.write("\t")
            f.write("\n")


def remove_redundant(target_seqs):
    seq_ids = [record.id for record in target_seqs]
    dup_seqs = {k[0]: [] for k in [seq_name for seq_name in list(Counter(seq_ids).items()) if seq_name[1] > 1]}

    count = 0
    while count < len(target_seqs):
        try:
            dup_seqs[target_seqs[count].id].append(target_seqs[count])
            del target_seqs[count]
        except KeyError:
            count += 1

    for record in dup_seqs.keys():
        target_seqs.append(dup_seqs[record][0])


def get_indv_genes(seq_list):
    seq_dict = {}
    for record in seq_list:
        temp_gn = record.id.split(".")[-1]
        seq_dict.setdefault(temp_gn, [])

        seq_dict[temp_gn].append(record)
    return seq_dict


if __name__ == '__main__':
    seq_records = list(SeqIO.parse(args.scraped_genes_fasta, "fasta"))
    remove_redundant(seq_records)

    # remove excluded loci from the input file
    if args.exclude != "":
        exclude_ids = [x[0] for x in utils.parse_file_nohead_tolist(args.exclude)]

        x = 0
        while x < len(seq_records):
            if seq_records[x].id in exclude_ids:
                print(seq_records[x])
                del seq_records[x]
                x -= 1
            x += 1

    tax_names = set(k.description.replace(k.id, "").replace(k.id.split(".")[-1], "").strip() for k in seq_records)
    if args.pre_fmt:
        tax_names = set(k.split("DESCRIPTION")[0].strip() for k in tax_names)

    uniq_spp = sorted(list(tax_names))

    aberrant_spp = clean_species(uniq_spp)
    novels, subs = get_weird_spp(uniq_spp)
    # uniq_spp = [x for x in uniq_spp if x not in novels and x not in subs]

    uniq_spp = [x for x in uniq_spp if x not in novels]
    # filter out undefined species
    uniq_spp = [x for x in uniq_spp if " sp." not in x]

    novels_translation = {x: "" for x in novels}
    novel_seqs = []
    if not args.pre_fmt:
        # reformat novel species so supercrunch can understand them and update sequence records
        for sp in novels:
            # get the specific name from the whole name
            new_sp_name = "_".join(sp.split(" ")[1:]).replace(".", "")
            # add the generic name back to the specific
            new_sp_name = " ".join([sp.split(" ")[0], new_sp_name])
            # replace any oddball characters with normal characters
            new_sp_name = new_sp_name.replace(":", "-").replace("/", "-").replace("\'", "").replace("\"", "")

            # update translation dictionary
            novels_translation[sp] = new_sp_name

        # update sequence dictionary
        k = 0
        while k < len(seq_records):
            temp_gene = seq_records[k].description.split(" ")[-1]
            temp_desc = " ".join(seq_records[k].description.split(" ")[:-1])
            seq_records[k].description = "%s DESCRIPTION %s" % (temp_desc, temp_gene)

            temp_sp = seq_records[k].description.replace(seq_records[k].id, "")
            temp_sp = temp_sp.replace(seq_records[k].id.split(".")[-1], "").strip()
            temp_sp = temp_sp.replace("DESCRIPTION", "").strip()

            # reformat temp spp if they're to remain in the output
            if args.include_odds:
                # temp_sp = temp_sp.split("DESCRIPTION:")[0].strip()
                if temp_sp not in novels and temp_sp not in uniq_spp:
                    del seq_records[k]
                    continue
                if temp_sp in novels:
                    seq_records[k].description = seq_records[k].description.replace(temp_sp,
                                                                                    novels_translation[temp_sp])
                    temp_sp = novels_translation[temp_sp]
                # remove space between specific and generic
                temp_sp = temp_sp.replace(" ", "_")
            # remove temp spp if they're not supposed to remain
            else:
                if temp_sp not in uniq_spp:
                    print(seq_records[k])
                    print("-----------")
                    novel_seqs.append(seq_records[k])
                    del seq_records[k]
                    continue

            # update sequence records to contain species names
            """
            seq_records[k].id = "%s.%s" % (temp_sp, seq_records[k].id)
            seq_records[k].name = "%s.%s" % (temp_sp, seq_records[k].name)
            seq_records[k].description = "%s.%s" % (temp_sp, seq_records[k].description)
            """
            k += 1

    if args.collapse_subs:
        uniq_spp = [x for x in uniq_spp if x not in subs]
        sub_translate = {x: " ".join(x.split(" ")[:2]) for x in subs}
        spp_w_subs = sorted(list(set(sub_translate.values())))
        include_subspp = set()
        for sp in spp_w_subs:
            try:
                nominate_sub = "%s %s" % (sp, sp.split(" ")[-1])
                sub_translate[nominate_sub]
                include_subspp.add(nominate_sub)
            except KeyError:
                for subsp in sub_translate.keys():
                    if subsp.startswith(sp):
                        include_subspp.add(subsp)
        # remove the subspecies that aren't nominate or don't represent the species
        pprint(sub_translate)
        pprint(include_subspp)
        for x in set(sub_translate.keys())-include_subspp:
            del sub_translate[x]

        uniq_spp.extend(sub_translate.values())
        uniq_spp = sorted(list(set(uniq_spp)))
        record = 0
        while record < len(seq_records):
            tmp_spname = seq_records[record].description.replace(seq_records[record].id, "")
            tmp_spname = tmp_spname.replace(seq_records[record].id.split(".")[-1], "").strip()
            tmp_spname = tmp_spname.split("DESCRIPTION")[0].strip()

            if tmp_spname in include_subspp:
                seq_records[record].description = seq_records[record].description.replace(tmp_spname,
                                                                                          sub_translate[tmp_spname])
            # delete sequence records for non-nominate or representative subspecies
            elif len(tmp_spname.split(" ")) == 3 and (tmp_spname not in include_subspp):
                print(seq_records[record])
                print("-----------")
                del seq_records[record]
                record -= 1
            record += 1

    # remove abberant spp from the fasta file


    gene_dict = get_indv_genes(seq_records)
    fasta_dir = "Parsed-Fasta-Files"
    gene_counts = []

    try:
        os.mkdir(fasta_dir)

    except FileExistsError:
        pass

    for gene_name in gene_dict.keys():
        if gene_dict[gene_name] != []:
            # print(gene_dict[gene_name])
            SeqIO.write(gene_dict[gene_name], "%s/%s.fasta" % (fasta_dir, gene_name), "fasta")
        gene_counts.append([gene_name, len(gene_dict[gene_name])])

    with open("gene_counts.txt", "w") as count_file:
        for line in gene_counts:
            for ele in line:
                count_file.write(str(ele))
                count_file.write("\t")
            count_file.write("\n")

    # write taxonomy files
    SeqIO.write(seq_records, "%s_tax-reformat.fasta" % Path(args.scraped_genes_fasta).stem, "fasta")

    with open("%s-%s" % (dv.SPP_FILE_NAME, "normal_taxa"), 'w') as species_file:
        for x in uniq_spp:
            species_file.write(x)
            species_file.write("\n")

    with open("%s-%s" % (dv.SPP_FILE_NAME, "odd_taxa"), 'w') as species_file:
        if not args.pre_fmt:
            for x in novels_translation.keys():
                species_file.write(x)
                species_file.write("\t")
                species_file.write(novels_translation[x])
                species_file.write("\n")

    if not args.pre_fmt:
        with open("%s-%s" % (dv.SPP_FILE_NAME, "combined_taxa"), 'w') as species_file:
            for x in uniq_spp:
                species_file.write(x)
                species_file.write("\n")
            for x in novels_translation.keys():
                species_file.write(novels_translation[x])
                species_file.write("\n")

    # nov_spp_dict = make_gene_dict(novels, uniq_genes, g_list)
    # sub_spp_dict = make_gene_dict(subs, uniq_genes, g_list)
