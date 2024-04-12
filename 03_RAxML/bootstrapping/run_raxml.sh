input_align=../pyron_final_alignment.phy
input_part=../pyron_8_partition.txt
run_name=pyron_final_rnd3

# Bootstrapping analysis
output_dir=/home/astewart/pyron_final_raxml/rnd3_100bs/
mkdir -p ${output_dir}
raxmlHPC-PTHREADS-SSE3 -s ${input_align} -q ${input_part} -n ${run_name}_BS.tre -m GTRCAT -T 31 -x $RANDOM -N 100 -p $RANDOM -k -w ${output_dir} &> ${run_name}_BS.log &
