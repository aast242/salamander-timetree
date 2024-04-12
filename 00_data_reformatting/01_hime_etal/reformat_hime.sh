# running the following command should make a file called 'hime_ahe_reformatted.fasta'
# that will be identical to 'hime_ahe_reformatted.fasta' in the starting files for SuperCRUNCH
hime_ahe_reformat.py ./inputs/KEBL01.1.gbff ./inputs/focal_hime_spp --add_dict ./inputs/hime-rovito-shen_gene_dict_scraped_sequences.txt
