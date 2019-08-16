#!/bin/bash
#SBATCH --job-name=picard_sort
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
java -Xmx8G -jar $PICARD SortSam \
	INPUT=na12878.sam \
	OUTPUT=na12878_sort.sam \
	SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=SILENT 

date

