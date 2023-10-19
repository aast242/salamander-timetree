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

parser = argparse.ArgumentParser(description="Program: Hime et al. AHE Reformatting\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s gene_file [options]')
parser.add_argument("ahe_gbff", help='A GenBank formatted file containing the sequences'
                                     ' from Hime et al. (2021)', type=str)
# this file can be downloaded at the link below:
# https://www.ncbi.nlm.nih.gov/Traces/wgs/?display=download
parser.add_argument("focal_taxa", help='File containing the species to pull from the Hime et al. data', type=str)


args = parser.parse_args()


if __name__ == '__main__':
    raw_hime_data = SeqIO.parse(args.ahe_gbff, "genbank")

    focal_species = [k[0] for k in utils.parse_file_nohead_tolist(args.focal_taxa)]

    """
    salamander_spp = set()
    for record in raw_hime_data:
        if "Caudata" in record.annotations["taxonomy"]:
            salamander_spp.add(record.annotations["organism"])
    salamander_spp = sorted(list(salamander_spp))
    with open("hime_salamander_spp.txt", 'w') as f:
        for sp in salamander_spp:
            f.write(sp)
            f.write("\n")
    """

    hime_rois = []
    for record in raw_hime_data:
        if record.annotations["organism"] in focal_species:
            seq_name = record.description.replace("TLS: %s " % record.annotations["organism"], "")
            seq_name = " ".join(seq_name.split(";")[0].split(" ")[1:])
            seq_name = seq_name.replace("anchored hybrid enrichment locus ", "HLL_")

            seq_id = "%s.%s" % (record.id, seq_name)

            temp_record = SeqRecord.SeqRecord(record.seq,
                                              description="%s %s" % (record.annotations["organism"], seq_name),
                                              name="%s %s" % (record.annotations["organism"], seq_name),
                                              id=seq_id)
            hime_rois.append(temp_record)
    SeqIO.write(hime_rois, "hime_ahe_reformatted.fasta", "fasta")
