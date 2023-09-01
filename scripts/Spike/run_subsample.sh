#!/bin/bash

# params
bampath=$1
samtools_path=$2
output_dirpath=$3

threads=10

prefix=`basename ${bampath} .bam`

ratio=$(cat ${output_dirpath}/${prefix}.ratio.txt)
${samtools_path} view -h -@ ${threads} -s ${ratio} -b -o ${output_dirpath}/${prefix}"_"${ratio}.bam ${bampath}
${samtools_path} index ${output_dirpath}/${prefix}"_"${ratio}.bam
