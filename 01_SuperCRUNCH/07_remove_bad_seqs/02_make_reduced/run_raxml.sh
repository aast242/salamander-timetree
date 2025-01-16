start_dir="$PWD"
mkdir -p gene_trees
for i in fastas/*.fasta; do	
	gene_name=$(basename ${i})
       	gene_name="${gene_name%%_oneseq*}";
	fasta_path=$(readlink -f "${i}")
	cd gene_trees
	echo $gene_name
	echo $fasta_path
	echo "---------"
	mkdir -p ${gene_name}
	cd $gene_name
	raxmlHPC-PTHREADS-SSE3 -s ${fasta_path} -n ${gene_name}_raxml.out -m GTRCAT -f c -T 1 -x $RANDOM -N 100 -p $RANDOM -k &
	cd $start_dir
done
