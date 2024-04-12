# running the following command will create a directory called "MrBayes_nexus_Files_SC-reformat" with
# the original nexus files reformatted as fasta files and ready for SuperCRUNCH.
# catting all of the files together (done by the setup script) 
# will make a file called "williams_2013_markers_refmt.fasta" that is identical 
# to "williams_2013_markers_refmt.fasta" in the starting SuperCRUNCH directory
williams_reformat.py MrBayes_nexus_Files/ taxonomy_transtable.txt
