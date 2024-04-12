mkdir -p blastn_out6
for i in ../outgroup_genic/fastas/*genic_regions.fasta; do 
	db_base=$(basename ${i});
	echo $db_base;
	blastn -subject $i -query pyron_williams_hime_genbank.fasta -outfmt 6 -task dc-megablast \
		-out ./blastn_out6/${db_base}.out6 -num_threads 4 -evalue 1e-10; 
done
