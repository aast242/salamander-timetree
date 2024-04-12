# to run on the combined scraped gbk sequences from the search performed in the paper
# the 'combined-scraped_genes.fasta' file output by this command will be identical to 
# 'hime-rovito-shen-scraped_genes.fasta' in the starting files for the SuperCRUNCH pipeline
gene_scraper_new.py ./combined.gbk --skip_seqpull --gbk_file \
	--add_dict ./inputs/hime-rovito-shen_gene_dict_scraped_sequences.txt \
	--taxid ./inputs/amphibian_outgroup_taxids.txt --full_mito

# If you'd like to re-do the GenBank scrape, uncomment the following lines!
#gene_scraper_new.py ./inputs/hime-rovito-shen.txt \
#        --add_dict ./inputs/hime-rovito-shen_gene_dict_scraped_sequences.txt \
#        --taxid ./inputs/amphibian_outgroup_taxids.txt --full_mito
