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
parser.add_argument("gene_file", help='File containing genes to search for', type=str)
parser.add_argument("--taxid", default="8293", help="the txid for the taxon you're pulling genes from", type=str)
parser.add_argument("--full_mito", action="store_true", help="performs a search for full mitochondrial genomes")
parser.add_argument("--include_mrna", action="store_true", help="allows mRNAs in the search results")
parser.add_argument("--print_products", action="store_true",
                    help="ONLY print the unique products in the gbk file")
parser.add_argument("--gbk_file", action="store_true",
                    help="interpret the gene_file as a gbk file and start later in the process")
parser.add_argument("--gbk_dir", type=str, default="",
                    help="interpret the gene_file as a directory containing gbk files and write a file that can "
                         "be used as a dictionary")
parser.add_argument("--skip_seqpull", action="store_true", help="skips downloading sequences from GenBank")
parser.add_argument("--add_dict", type=str, default="",
                    help="a file containing translations for gene names made using the \'--gbk_dir\' option"
                         " and !!CHECKED FOR INCONSISTENCIES!!")

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delimited on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

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


def create_queries(file_dict):

    if Path(args.taxid).exists():
        taxon_list = [i[0] for i in utils.parse_file_nohead_tolist(args.taxid)]
    else:
        taxon_list = [args.taxid]

    constructed_queries = {}
    extra_queries = []
    for goi in file_dict.keys():
        constructed_queries.setdefault(goi, [])
        for taxon in taxon_list:
            temps = []
            for ele in range(0, len(file_dict[goi])):
                if file_dict[goi][ele] != "":
                    if " " in file_dict[goi][ele]:
                        temps.append("(" + "[Title] AND ".join(file_dict[goi][ele].split(" ")) + "[Title])")
                    else:
                        temps.append("(%s[Title] OR %s[Gene Name])" % (file_dict[goi][ele], file_dict[goi][ele]))

            # combine the individual pieces and add the taxon requirement
            #finals = "(" + " OR ".join(temps) + ") AND txid%s[organism:exp] AND ddbj_embl_genbank[filter]" % args.taxid
            finals = "(" + " OR ".join(temps) + ") AND txid%s[organism:exp]" % taxon
            if not args.include_mrna:
                finals += " AND biomol_genomic[PROP]"
            constructed_queries[goi].append(finals)

    if args.full_mito:
        constructed_queries.setdefault("mitochondria", [])
        for taxon in taxon_list:
            constructed_queries["mitochondria"].append("((mitochondrial OR mitochondrion OR mitochondria) "
                                                       "AND genome[Title]) AND txid%s[organism:exp] "
                                                       "AND biomol_genomic[PROP] " % taxon)

    #print(constructed_queries)
    #exit()

    return constructed_queries


def dump_acc_to_gbk(query_list):
    loop_len = len(query_list)
    loop_iter = 1
    for goi in query_list.keys():
        print("Starting GOI %s..." % goi)
        if Path("%s/%s.gbk" % (dv.OUTDIR_PREFIX, goi)).exists():
            print("File for %s already exists!" % goi)
            print("Skipping...\n")
            loop_iter += 1
            continue
        print("Searching nucleotide database with each taxid. Here is an example search: ")
        print(query_list[goi][0])
        compiled_results = []
        for taxon in query_list[goi]:
            # print("Search command: %s" % taxon)
            # print("Finding how many search results...")
            # get the number of results returned by the search
            handle = Entrez.esearch(db="nucleotide", term=taxon)
            record = Entrez.read(handle)
            count = record['Count']
            handle.close()

            # print("%s search results!" % count)
            # print("Fetching them from GenBank...")

            # actually get everything searched
            handle = Entrez.esearch(db="nucleotide", retmax=count, term=taxon)
            # handle = Entrez.esearch(db="nucleotide", retmax=5, term=taxon)
            record = Entrez.read(handle)
            record = list(record["IdList"])
            compiled_results.extend(record)
            handle.close()
        print("Finished searches for GOI %s!" % goi)
        if len(compiled_results) == 1:
            print("There is %s retrieved record for %s" % (len(compiled_results), goi))
            print("Pulling GenBank record...")
        else:
            print("There are %s retrieved records for %s" % (len(compiled_results), goi))
            print("Pulling GenBank records...")
        record = list(utils.chunks(compiled_results, 500))

        chunk_count = 0
        for chunk in record:
            handle = Entrez.efetch(db="nucleotide", id=chunk, rettype="gb", retmode="gbwithparts")

            with open("%s/chunk_%s.gbk" % (dv.OUTDIR_PREFIX, chunk_count), "w") as gbk_file:
                for i in list(handle):
                    gbk_file.write(str(i))
            chunk_count += 1

            if chunk_count % 2 == 0:
                if chunk_count * 500 > len(compiled_results):
                    print("processed all samples")
                else:
                    print("processed %s samples" % (chunk_count * 500))
                # temp_prods = [i for i in temp_prods if "cytochrome c oxidase" in i.lower()]
        # combine gbk chunks into one file
        temp_gbk_chunks = ["%s/chunk_%s.gbk" % (dv.OUTDIR_PREFIX, i) for i in range(0, chunk_count)]

        with open("%s/%s.gbk" % (dv.OUTDIR_PREFIX, goi), 'wb') as final_file:
            for gbk in temp_gbk_chunks:
                with open(gbk, 'rb') as temp_file:
                    shutil.copyfileobj(temp_file, final_file)
                os.remove(gbk)
        print("finished with goi (%s/%s)\n" % (loop_iter, loop_len))
        loop_iter += 1
    print("finished with gene products!")


def parse_gbk_entry(gbk_record, genes_of_interest):
    # TODO: based on these analyses, there are a TON of "misc features" that have 12s and 16s fragments
    # maybe eventually try to pull those out? probably not worth it, but could be interesting
    # print(gbk_record.features)
    feature_list = []
    for feature in gbk_record.features:
        try:
            product_identity = feature.qualifiers["product"][0].lower()
            try:
                # Make sure this is a gene we actually care about!
                feature_gene = genes_of_interest[product_identity]

                feature_acc = gbk_record.id
                feature_organism = gbk_record.annotations["organism"]
                feature_range = feature.location
                feature_range = "%s:%s" % (int(feature_range.start), int(feature_range.end))
                feature_date = gbk_record.annotations["date"]
                feature_seq = feature.extract(gbk_record.seq)
                feature_length = len(feature_seq)

                feature_list.append([feature_organism, feature_acc, feature_range, feature_gene,
                                     feature_date, feature_length])

            except KeyError:
                pass
        except KeyError:
            pass
    return feature_list


def get_products_in_gbk(gbk_record):
    short_gn = ""
    long_gn = ""
    for feature in gbk_record.features:
        # products are only present for mRNAs and CDS
        try:
            product_identity = feature.qualifiers["product"][0].lower()
            long_gn = product_identity
        except KeyError:
            pass
        # genes are present for gene, mRNA, and CDS
        try:
            product_identity = feature.qualifiers["gene"][0].lower()
            short_gn = product_identity
        except KeyError:
            pass
    products = tuple([short_gn, long_gn])
    #if long_gn == "kiaa":
    #    print(gbk_record)
    return products


def make_gbk_dict(parsed_gbk):
    gbk_dict = {}
    for record in parsed_gbk:
        gbk_dict.setdefault(record.id, record)
    return gbk_dict


def merge_tuples(list_of_tuples):
    new_tuples = list_of_tuples.copy()
    new_tuples = [set(i) - {''} for i in new_tuples]
    merged_list = []
    i = 0
    while i < len(new_tuples):
        n = 0
        while n < len(new_tuples):
            if n == i:
                n += 1
                continue

            if new_tuples[i].intersection(new_tuples[n]) != set():
                new_tuples[i].update(new_tuples[n])
                del new_tuples[n]
                n -= 1
            n += 1
        merged_list.append(new_tuples[i])
        del new_tuples[i]
    merged_list = [sorted(list(i), key=len) for i in merged_list]
    return merged_list


if __name__ == '__main__':
    entry_file_name = Path(args.gene_file).stem
    if not (args.gbk_file or (args.gbk_dir != "")):
        if not Path(dv.OUTDIR_PREFIX).exists():
            os.mkdir(dv.OUTDIR_PREFIX)

        g_path = Path(args.gene_file).resolve()
        gene_names = utils.parse_file_nohead_todict(g_path, 0)
        expected_genes = create_queries(gene_names)
        if not args.skip_seqpull:
            dump_acc_to_gbk(expected_genes)

        combined_gbk_filename = "%s.gbk" % dv.COMB_FILE

        downloaded_genes = list(gene_names.keys())
        if args.full_mito:
            downloaded_genes.append("mitochondria")

        with open(combined_gbk_filename, 'wb') as comb_gbk:
            for gn in ["%s/%s.gbk" % (dv.OUTDIR_PREFIX, k) for k in downloaded_genes]:
                with open(gn, 'rb') as gene_gbk:
                    shutil.copyfileobj(gene_gbk, comb_gbk)

        print("Parsing combined GenBank file for all queried genes...")
        comb_gbk_file = SeqIO.parse(combined_gbk_filename, "genbank")
        gbk_dictionary = make_gbk_dict(comb_gbk_file)
    else:
        if args.gbk_dir != "":
            if not Path(args.gbk_dir).is_dir():
                print("FATAL: path provided to \'--gbk_dir\' is not a directory!")
                exit()

            orig_gene_file = utils.parse_file_nohead_tolist(args.gene_file)
            orig_gene_file = [[j.lower() for j in k if j != ""] for k in orig_gene_file]

            gbk_files = ["%s/%s" % (Path(args.gbk_dir).resolve(), k)
                         for k in os.listdir(Path(args.gbk_dir).resolve()) if k.endswith(".gbk")]
            gbk_dicts = {}
            for k in gbk_files:
                gene_name = Path(k).stem
                temp_gbk_dict = make_gbk_dict(SeqIO.parse(k, "genbank"))
                temp_prods = set()
                temp_prods.add(gene_name.lower())
                for acc in temp_gbk_dict.keys():
                    for j in get_products_in_gbk(temp_gbk_dict[acc]):
                        temp_prods.add(j.lower())
                temp_prods = sorted([j for j in list(temp_prods) if j], key=len)

                if temp_gbk_dict != {}:
                    gbk_dicts.setdefault(gene_name, temp_prods)
            for k in orig_gene_file:
                try:
                    gbk_dicts[k[0]]
                except KeyError:
                    gbk_dicts[k[0]] = sorted([j for j in k], key=len)
            dict_file_name = "%s_gene_dict_%s.txt" % (entry_file_name, Path(args.gbk_dir).stem)
            print("Writing gene dictionary to \'%s\'" % dict_file_name)
            print("MAKE SURE to manually check it!")
            with open(dict_file_name, "w") as f:
                for acc in sorted(list(gbk_dicts.keys())):
                    if acc != "mitochondria":
                        for names in gbk_dicts[acc]:
                            f.write("%s\t%s\n" % (acc.lower(), names.lower()))
                        f.write("\n")
            exit()

        elif args.gbk_file:
            print("Parsing combined GenBank file for all queried genes...")
            comb_gbk_file = SeqIO.parse(args.gene_file, "genbank")
            gbk_dictionary = make_gbk_dict(comb_gbk_file)

        if args.print_products:
            print("getting all unique product names in the GenBank file...")
            uniq_products = set()
            for acc in gbk_dictionary.keys():
                uniq_products.add(get_products_in_gbk(gbk_dictionary[acc]))
            uniq_products = sorted(list(uniq_products))
            uniq_products = merge_tuples(uniq_products)
            pprint(uniq_products)

            # prod_dict = {}
            # for goi in uniq_products:
            #     for k in range(0, len(goi)):
            #         prod_dict[goi[k]] = goi[0]
            # pprint(prod_dict, sort_dicts=False)
            exit()
            with open("%s-uniq_products.txt" % entry_file_name, 'w') as f:
                for line in uniq_products:
                    for entry in line:
                        f.write(entry)
                        f.write("\t")
                    f.write("\n")
            exit()

    print("Getting gene products from the parsed GenBank file...")
    all_features = []
    for acc in list(gbk_dictionary.keys()):
        temp_feats = parse_gbk_entry(gbk_dictionary[acc], dv.MITO_GOIS)
        all_features.extend(temp_feats)

    print("removing duplicate entries...")
    all_features = list(all_features for all_features, _ in itertools.groupby(all_features))

    print("checking for gapped sequences within accessions...")
    utils.check_gaps(all_features)
    all_features.sort(key=lambda k: (k[dv.ORDER_DICT["species"]], k[dv.ORDER_DICT["accession"]]))

    print("Pulling sequences...")
    all_seqrecords = []
    for line in all_features:
        temp_gbk_entry = gbk_dictionary[line[dv.ORDER_DICT["accession"]]]

        foc_range = [int(k) for k in line[dv.ORDER_DICT["range"]].split(":")]
        if foc_range[0] < foc_range[1]:
            seq_pol = +1
        else:
            seq_pol = -1
        foc_range = SeqFeature.FeatureLocation(min(foc_range), max(foc_range), strand=seq_pol)

        """
        long_desc = "label not included in dictionary"
        try:
            long_desc = GENE_DESCRIPTION[line[dv.ORDER_DICT["gene"]]]
        except KeyError:
            pass
        """

        all_seqrecords.append(
            SeqIO.SeqRecord(foc_range.extract(temp_gbk_entry.seq),
                            id="%s.%s" % (line[dv.ORDER_DICT["accession"]], line[dv.ORDER_DICT["gene"]]),
                            name="%s.%s" % (line[dv.ORDER_DICT["accession"]], line[dv.ORDER_DICT["gene"]]),
                            description="%s %s" % (line[dv.ORDER_DICT["species"]], line[dv.ORDER_DICT["gene"]],)
                            )
        )

    print("Writing files...")
    # rewrite the pared genbank table and remove the old table
    if not args.gbk_file:
        SeqIO.write([gbk_dictionary[k] for k in gbk_dictionary.keys()],
                    "%s-scraped_genes_table.gbk" % entry_file_name,
                    "genbank")
        os.remove(combined_gbk_filename)

    # write a fasta with the genic features
    SeqIO.write(all_seqrecords,
                "%s-scraped_genes.fasta" % entry_file_name,
                "fasta")

    # write the extracted genic features
    with open("%s-scraped_genes_table.txt" % entry_file_name, 'w') as f:
        for line in all_features:
            for entry in line:
                f.write(str(entry))
                f.write("\t")
            f.write("\n")
