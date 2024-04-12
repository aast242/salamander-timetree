start_dir="$PWD"
mkdir -p ./tree_imgs
cd ./tree_imgs

for i in $start_dir/gene_trees/*/RAxML_bipartitions.*; do
	python3 render_tree.py $i
done
