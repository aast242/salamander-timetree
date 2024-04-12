while IFS=$'\t' read -r -a rename_array
do
 old_name="${rename_array[0]}"
 new_name="${rename_array[1]}"
 sed -i "s/${old_name}/${new_name}/g" mito_simple_extracted.fasta 
done < rename_file 

make_taxfile.py mito_simple_extracted.fasta --collapse_subs --pre_fmt --exclude paramesotriton_exclude  

