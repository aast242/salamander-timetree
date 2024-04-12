#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")

python3 $SC_SCRIPTS/Align.py -i ${in_dir} -o ${out_dir} -a macse --mpath $SC_SCRIPTS/macse_v2.05.jar --table vertmtdna --pass_fail --mem 10
