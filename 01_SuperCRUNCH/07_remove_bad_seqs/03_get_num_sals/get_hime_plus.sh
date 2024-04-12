#!/bin/bash

in_dir=$(readlink -f "${1}")

for i in $in_dir/*.fasta; do
	og_count=$(grep -f focal_hime_spp -c $i)
	salamander_count=$(grep -i -f salamander_names -c $i)
	hime_count=$(grep "KEBL" -c $i)
	tot_count=$(grep ">" -c $i)
	og_plus_hime=$((${og_count} + ${hime_count}))

	if  [[ ${og_plus_hime} -lt ${tot_count} ]] 
	then
		true
		#echo $i >> hime_plus_fastas.txt
		#echo "Not just hime!"
	fi
	
	if  [[ ${salamander_count} -ge 10 ]]
	then
		printf "%s\t%s\t%s\t%s\t%s\n" $i $og_count $hime_count $tot_count $salamander_count >> enough_salamanders.txt
		#printf "Gene Name: %s\n" $i
		#printf "Outgroup Species Count: %s\n" $og_count
		#printf "Number of Hime species: %s\n" $hime_count
		#printf "Total Number of Sequences: %s\n" $tot_count
		#printf "Salamander Count: %s\n" $salamander_count
		#echo "------------"
	else		
		printf "%s\t%s\t%s\t%s\t%s\n" $i $og_count $hime_count $tot_count $salamander_count >> not_enough_salamanders.txt
		echo "~~~~~~NOT ENOUGH SALAMANDERS~~~~~~~"
		printf "Gene Name: %s\n" $i
		printf "Total Number of Sequences: %s\n" $tot_count
		printf "Salamander Count: %s\n" $salamander_count
		echo "------------"
	fi
	
done

