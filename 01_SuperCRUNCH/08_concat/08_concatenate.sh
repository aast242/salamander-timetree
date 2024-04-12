#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")
python3 $SC_SCRIPTS/Concatenation.py -i ${in_dir} -o ${out_dir} --informat fasta --outformat phylip --symbol dash
