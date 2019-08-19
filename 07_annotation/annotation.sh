#!/bin/bash
#SBATCH --job-name=annotation
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

module load samtools/1.7
bgzip ../06_variants/na12878_q20.recode.vcf
tabix ../06_variants/na12878_q20.recode.vcf.gz 


module unload  samtools/1.7


module load bcftools/1.6
bcftools annotate -c ID \
	-a ../02_reference_data/reference_chr20.vcf.gz ../06_variants/na12878_q20.recode.vcf.gz \
	> na12878_q20_annot.vcf


