#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alexander Stewart"
__copyright__ = "Copyright 2023"
__credits__ = ["Alexander Stewart"]
__license__ = "GPL3"
__maintainer__ = "Alexander Stewart"
__status__ = "Development"

# base python imports #
import argparse
import os
import itertools
import sys
from shutil import rmtree
from pprint import pprint
from pathlib import Path
from re import findall
from collections import Counter
from datetime import datetime

# external package imports #
from ete3 import PhyloTree
from Bio import Entrez, SeqIO, SeqRecord

# local imports #
from defaults import ProgDefaults as dv
import utils

parser = argparse.ArgumentParser(description="Program: Render Tree\n"
                                             "Version: 1.0\n",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s tree_file [options]')
parser.add_argument("tree_file", help="", type=str)

args = parser.parse_args()

outgroup_species = ["Latimeria_chalumnae", "Homo_sapiens", "Anolis_carolinensis", "Chrysemys_picta", "Gallus_gallus",
                    "Rhinatrema_bivittatum", "Microcaecilia_unicolor", "Geotrypetes_seraphini",
                    "Xenopus_tropicalis"]


tree_obj = utils.parse_file_nohead_tolist(args.tree_file)[0][0]
tree_obj = PhyloTree(tree_obj)

leaf_labels = [str(i).replace("\n--", "") for i in tree_obj.get_leaves()]
species_in_tree = ["_".join(i.split("_")[:2]) for i in leaf_labels]

curr_root = ""
for i in outgroup_species:
    if i in species_in_tree:
        for j in leaf_labels:
            if i in j:
                curr_root = j
                break
    if curr_root != "":
        break

if curr_root != "":
    tree_obj.set_outgroup(curr_root)

tree_obj.ladderize()
tree_obj.render("%s.svg" % Path(args.tree_file).stem)
tree_obj.render("%s.png" % Path(args.tree_file).stem)

