#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")

python3 $SC_SCRIPTS/Align.py -i ${in_dir} -o ${out_dir} -a mafft --threads 4 --accurate
