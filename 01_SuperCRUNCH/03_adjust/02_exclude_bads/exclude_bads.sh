for i in ./inputs/*.fasta; do
	echo $i
	remove_specified_seqs.py ${i} exclude_ids
done
