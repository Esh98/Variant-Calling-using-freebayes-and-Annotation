#!/bin/bash
#SBATCH --job-name=freebayes_variants
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

module load freebayes/1.1.0

freebayes -f ../02_reference_data/chr20.fa \
	../05_align/na12878_sort_marked.bam > na12878.vcf 


