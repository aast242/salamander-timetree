# running the following command will create a directory called "haplotypes_filtered_SC-reformat" with
# the original fasta files reformatted and ready for SuperCRUNCH.
# catting all of the files together (done by the setup script) 
# will make a file called "everson_2021_markers_refmt.fasta" that is identical 
# to "everson_2021_markers_refmt.fasta" in the starting SuperCRUNCH directory
everson_reformat.py ./haplotypes_filtered/ new_taxonomy_transtable.txt
