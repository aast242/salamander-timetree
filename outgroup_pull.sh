#!/usr/bin/env bash

tsv_dir=$1
fileend=.tsv
script_dir=$(dirname "$0")
pull_prog=outgroup_pull_genes.py
comb_file=combined_outgroups.genes
dict_add=hime-rovito-shen_gene_dict_scraped_sequences.txt

for i in ${tsv_dir}*${fileend}; do
       	python3 ${script_dir}/${pull_prog} ${i} --silent --add_dict ${script_dir}/${dict_add}
	if [ $? -ne 0 ]; then
		echo "Something went wrong! Try running the script again."
		break
	fi
done

# remove the combined file if it's there
if [ -e ./${comb_file} ]; then
	rm ./${comb_file}
fi

# cat
gene_files=(./*.genes)
head -n1 ${gene_files[0]} >> ${comb_file}

for i in "${gene_files[@]}"; do
	tail -n+2 $i >> ${comb_file}
done

