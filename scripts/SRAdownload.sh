#!/bin/bash

VAR=$(cat $1)

for i in ${VAR}
do
	echo "Downloading SRA entry:${i}"
    
	fastq-dump --gzip --defline-qual '+' ${i}

	echo "Done downloading ${i}"
done

