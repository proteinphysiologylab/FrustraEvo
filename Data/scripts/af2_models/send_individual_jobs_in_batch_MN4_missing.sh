#!/bin/bash
fasta_dir="$1"
output_models="$2"


for i in $(ls ${fasta_dir}/*.fa)
do
j=$(basename $i .fa)
output_dir=${output_models}/${j}

if [ ! -f ${output_dir}/*_relaxed.pdb ];then
    if cat $i | grep -v '>' | grep -q 'X'; then
	    echo $i 'contains_X'	
    else
    	    echo $i 'unfinished model'
	    sbatch run_AF_minotauro_support_MN4.sh $i ${output_dir}
    fi
fi
done

