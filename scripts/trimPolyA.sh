#!/bin/bash

for i in *trimmed.fq.gz
do
    java -jar /usr/share/java/trimmomatic-0.39.jar SE ${i} trimmed/`basename ${i} .fq.gz`_polyA.fq.gz ILLUMINACLIP:/usr/share/trimmomatic/PolyA.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	

    fastqc -o trimmed/fastqc_final/ trimmed/`basename ${i} .fq.gz`_polyA.fq.gz
done

