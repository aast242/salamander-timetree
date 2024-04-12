#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")

python3 $SC_SCRIPTS/Cluster_Blast_Extract.py -i ${in_dir} -o ${out_dir} -b dc-megablast --threads 4
