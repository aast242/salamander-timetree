#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")
python3 $SC_SCRIPTS/Fasta_Filter_by_Min_Seqs.py -i ${in_dir} -o ${out_dir} --min_seqs 15


in_dir=$(readlink -f "./num_seq_adjust/outputs/Filtered-Fasta-Files/")
out_dir=$(readlink -f "./direction_adjust/")
#python3 $SC_SCRIPTS/Adjust_Direction.py -i ${in_dir} -o ${out_dir} --accurate --threads 6


in_dir=$(readlink -f "./translation_adjust/simple/inputs/")
out_dir=$(readlink -f "./translation_adjust/simple/outputs/")
#python3 $SC_SCRIPTS/Coding_Translation_Tests.py -i ${in_dir} -o ${out_dir} --table standard --quiet

in_dir=$(readlink -f "./translation_adjust/mito/inputs/")
out_dir=$(readlink -f "./translation_adjust/mito/outputs/")
#python3 $SC_SCRIPTS/Coding_Translation_Tests.py -i ${in_dir} -o ${out_dir} --table vertmtdna --quiet
