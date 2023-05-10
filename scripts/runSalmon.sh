#!/bin/bash

GEO=$(cat geo_accessions.txt)


for i in ${GEO}
do
	SRR=$(grep ${i} SraRunTable.txt | cut -d ',' -f 1)
	SRR=${SRR}_trimmed_polyA.fq.gz
	
	salmon quant -i gencode_index -l A -o ${i} -r ${SRR}	
done

	

