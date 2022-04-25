#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-122

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_files.txt)

module load fastqc

fastqc --nogroup --outdir ./ $INPUT_FILE
