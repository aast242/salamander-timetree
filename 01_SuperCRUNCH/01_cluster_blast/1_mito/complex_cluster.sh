#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")
multi_dir=$(readlink -f ".")

python3 $SC_SCRIPTS/Reference_Blast_Extract.py -i ${in_dir} -o ${out_dir} --multisearch ${multi_dir}/multifile -b dc-megablast -m span --threads 8
