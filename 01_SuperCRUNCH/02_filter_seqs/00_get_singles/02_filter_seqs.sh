#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")
taxon_file=$(readlink -f "./captured_caudate_taxa-normal_taxa")

python3 $SC_SCRIPTS/Filter_Seqs_and_Species.py -i ${in_dir} -o ${out_dir} -s oneseq -f length -m 150 -t ${taxon_file} 
