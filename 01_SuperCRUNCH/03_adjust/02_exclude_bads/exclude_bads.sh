for i in ./inputs/*.fasta; do
	echo $i
	../../../dtop_saltree/salamander-timetree-main/remove_specified_seqs.py ${i} exclude_ids
done
