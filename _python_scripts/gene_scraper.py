#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

import argparse
from Bio import Entrez
from pathlib import Path
import itertools
import sys
from shutil import rmtree
from pprint import pprint
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
parser.add_argument("--txid", default="8293", help="the txid for the taxon you're pulling genes from", type=str)
parser.add_argument("--full_mito", action="store_true", help="performs a search for full mitochondrial genomes")
parser.add_argument("--include_mrna", action="store_true", help="allows mRNAs in the search results")

# Gene file is organized so that each gene is on a line and synonyms for that gene are tab delimited on the line. eg,
# COX1\tCO1\tcytochrome oxidase subunit 1
# COX2\tCO2\tcytochrome oxidase subunit 2

args = parser.parse_args()


def create_queries(file_dict):
    constructed_queries = []
    extra_queries = []
    for goi in file_dict.keys():
        temps = []
        # set est_seq_len large so nothing is excluded at first
        for ele in range(0, len(file_dict[goi])):
            if " " in file_dict[goi][ele]:
                temps.append("(" + "[Title] AND ".join(file_dict[goi][ele].split(" ")) + "[Title])")
            else:
                temps.append("(%s[Title] OR %s[Gene Name])" % (file_dict[goi][ele], file_dict[goi][ele]))

        # combine the individual pieces and add the taxon requirement
        #finals = "(" + " OR ".join(temps) + ") AND txid%s[organism:exp] AND ddbj_embl_genbank[filter]" % args.txid
        finals = "(" + " OR ".join(temps) + ") AND txid%s[organism:exp]" % args.txid
        if not args.include_mrna:
            finals += " AND biomol_genomic[PROP]"

        constructed_queries.append(finals)

    if args.full_mito:
        #constructed_queries.append("((mitochondrial OR mitochondrion OR mitochondria) AND genome[Title]) "
        #                           "AND txid%s[organism:exp] AND biomol_genomic[PROP] "
        #                           "AND ddbj_embl_genbank[filter]" % args.txid)
        constructed_queries.append("((mitochondrial OR mitochondrion OR mitochondria) AND genome[Title]) "
                                   "AND txid%s[organism:exp] AND biomol_genomic[PROP] " % args.txid)

    #print(constructed_queries)
    #exit()

    return constructed_queries


def get_genes(query_list):
    gene_products = []
    all_gene_names = set()
    for goi in query_list:
        print("Search command: %s" % goi)
        print("Finding how many search results...")
        # get the number of results returned by the search
        handle = Entrez.esearch(db="nucleotide", term=goi)
        record = Entrez.read(handle)
        count = record['Count']
        handle.close()

        print("%s search results!\nFetching them from GenBank..." % count)

        # actually get everything searched
        handle = Entrez.esearch(db="nucleotide", retmax=count, term=goi)
        #handle = Entrez.esearch(db="nucleotide", retmax=5, term=goi)
        record = Entrez.read(handle)
        record = list(utils.chunks(record["IdList"], 500))
        handle.close()
        chunk_count = 0
        for chunk in record:

            handle = Entrez.efetch(db="nucleotide", id=chunk,
                                   rettype="gb", retmode="xml")

            fetched = list(Entrez.read(handle))
            for sample in fetched:
                #if sample["GBSeq_organism"] == "Bolitoglossa splendida":
                #    print(sample)
                #    print(sample["GBSeq_feature-table"])
                for r in sample["GBSeq_feature-table"]:
                    try:
                        for feat in r["GBFeature_quals"]:
                            #if feat["GBQualifier_name"] == "product":
                            #    all_gene_names.add(feat["GBQualifier_value"].lower())
                            if feat["GBQualifier_name"] == "product" \
                                    and feat["GBQualifier_value"].lower() in dv.MITO_GOIS.keys():
                                gene_products.append([sample["GBSeq_organism"],
                                                      r["GBFeature_intervals"][0]["GBInterval_accession"],
                                                      "%s:%s" % (r["GBFeature_intervals"][0]["GBInterval_from"],
                                                                 r["GBFeature_intervals"][0]["GBInterval_to"]),
                                                      dv.MITO_GOIS[feat["GBQualifier_value"].lower()],
                                                      sample["GBSeq_create-date"]])
                                #if sample["GBSeq_organism"] == "Bolitoglossa splendida":
                                #    print(feat["GBQualifier_value"])
                                #    print(r)

                    except KeyError:
                        pass
            chunk_count += 1
            if chunk_count % 2 == 0:
                print("processed %s samples" % (chunk_count * 500) )
                # temp_prods = [i for i in temp_prods if "cytochrome c oxidase" in i.lower()]
        print("finished with goi\n")

    print("finished with gene products!")
    # pprint(sorted(list(all_gene_names)))

    return gene_products


if __name__ == '__main__':
    g_path = Path(args.gene_file).resolve()
    gene_names = utils.parse_file_nohead_todict(g_path, 0)
    expected_genes = create_queries(gene_names)
    organized_genes = get_genes(expected_genes)

    # add length to each entry
    for k in range(0, len(organized_genes)):
        len_pair = [int(j) for j in organized_genes[k][2].split(":")]
        total_len = max(len_pair) - min(len_pair)
        organized_genes[k].append(str(total_len))

    # remove duplicates
    # print(len(organized_genes))
    organized_genes.sort()
    organized_genes = list(organized_genes for organized_genes,_ in itertools.groupby(organized_genes))
    # print(len(organized_genes))

    with open("pulled_genes", 'w') as f:
        for line in organized_genes:
            for entry in line:
                f.write(entry)
                f.write('\t')
            f.write('\n')
