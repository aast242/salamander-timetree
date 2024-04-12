#!/bin/bash

in_dir=$(readlink -f "./inputs/")
out_dir=$(readlink -f "./outputs/")

python3 $SC_SCRIPTS/Fasta_Relabel_Seqs.py -i ${in_dir} -o ${out_dir} -r species
python3 $SC_SCRIPTS/Fasta_Relabel_Seqs.py -i ${in_dir} -o ${out_dir} -r species_acc
python3 $SC_SCRIPTS/Fasta_Relabel_Seqs.py -i ${in_dir} -o ${out_dir} -r accession
