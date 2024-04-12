for i in ./exclude_again/*.fasta; do
	echo $i
	remove_specified_seqs.py ${i} new_excludes
done
