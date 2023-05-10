#!/bin/bash

VAR=$(cat $1)

for i in ${VAR}
do
    echo "Downloading SRA entry:${i}"

    fastq-dump --gzip --defline-qual '+' ${i}

    echo "Done downloading ${i}"
done


for j in *fastq.gz
do
    java -jar /usr/share/java/trimmomatic-0.39.jar SE ${j} `basename ${j} .fastq.gz`_trimmed.fq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    fastqc -o fastqc_mid/ `basename ${j} .fastq.gz`_trimmed.fq.gz
done

for k in *trimmed.fq.gz
do
    java -jar /usr/share/java/trimmomatic-0.39.jar SE ${k} trimmed/`basename ${k} .fq.gz`_polyA.fq.gz ILLUMINACLIP:/usr/share/trimmomatic/PolyA.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


    fastqc -o trimmed/fastqc_final/ trimmed/`basename ${k} .fq.gz`_polyA.fq.gz
done

