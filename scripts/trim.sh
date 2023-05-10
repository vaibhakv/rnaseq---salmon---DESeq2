#!/bin/bash

for i in *fastq.gz
do
	java -jar /usr/share/java/trimmomatic-0.39.jar SE ${i} `basename ${i} .fastq.gz`_trimmed.fq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
 
	fastqc `basename ${i} .fastq.gz`_trimmed.fq.gz
done

