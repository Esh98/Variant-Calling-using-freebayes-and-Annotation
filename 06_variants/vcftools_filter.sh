#!/bin/bash
#SBATCH --job-name=vcftools_filter
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load vcftools/0.1.16

vcftools --vcf na12878.vcf \
	--minQ 20 \
	--recode \
	--recode-INFO-all \
	--out na12878_q20

