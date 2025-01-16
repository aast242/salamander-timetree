#!/bin/bash

#in_dir=$(readlink -f "./simple/inputs/")
#out_dir=$(readlink -f "./simple/outputs/")
#python3 $SC_SCRIPTS/Coding_Translation_Tests.py -i ${in_dir} -o ${out_dir} --table standard

in_dir=$(readlink -f "./mito/inputs/")
out_dir=$(readlink -f "./mito/outputs/")
python3 $SC_SCRIPTS/Coding_Translation_Tests.py -i ${in_dir} -o ${out_dir} --table vertmtdna
