#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

module load trimgalore

for file in *fastq.gz;do

        cutadapt -e 0.1 -q 20 -m 20 -O 1 -a NNAGATCGGAAGAGCACAC -a AGATCGGAAGAGCACAC -a ATCGGAAGAGCACAC -o ${file%.fastq.gz}.trimmed_reads.fastq.gz $file;

done
