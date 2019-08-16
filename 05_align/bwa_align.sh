#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bwa/0.7.17
bwa mem -M -t 2 \
	../02_reference_data/chr20 \
	../04_trimmed_reads/trimmed_NA12878_20_paired_1.fq ../04_trimmed_reads/trimmed_NA12878_20_paired_2.fq \
	-o na12878.sam 


