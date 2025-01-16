#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")
python3 $SC_SCRIPTS/Adjust_Direction.py -i ${in_dir} -o ${out_dir} --accurate --threads 12 
