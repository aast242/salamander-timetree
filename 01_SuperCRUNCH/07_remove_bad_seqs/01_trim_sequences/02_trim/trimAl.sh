#!/bin/bash

in_dir=$(readlink -f "../01_rename/outputs/Relabeled-by-species_acc/Relabeled-fasta-files")
out_dir=$(readlink -f "./species_acc")
python3 $SC_SCRIPTS/Trim_Alignments_Trimal.py -i ${in_dir} -o ${out_dir} -a both --gt 0.10 -f fasta

in_dir=$(readlink -f "../01_rename/outputs/Relabeled-by-species/Relabeled-fasta-files")
out_dir=$(readlink -f "./only_species")
python3 $SC_SCRIPTS/Trim_Alignments_Trimal.py -i ${in_dir} -o ${out_dir} -a both --gt 0.10 -f fasta

in_dir=$(readlink -f "../01_rename/outputs/Relabeled-by-accession/Relabeled-fasta-files")
out_dir=$(readlink -f "./accession")
python3 $SC_SCRIPTS/Trim_Alignments_Trimal.py -i ${in_dir} -o ${out_dir} -a both --gt 0.10 -f fasta
