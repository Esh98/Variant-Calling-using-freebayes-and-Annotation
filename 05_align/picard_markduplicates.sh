#!/bin/bash
#SBATCH --job-name=picard_markduplicates
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load picard/2.9.2

java -Xmx8G -jar $PICARD MarkDuplicates \
	INPUT=na12878_sort.sam \
	OUTPUT=na12878_sort_marked.bam \
	METRICS_FILE=metrics.txt \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT 

