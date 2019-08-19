#!/bin/bash
#SBATCH --job-name=snpEff
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=mcbstudent
#SBATCH --partition=mcbstudent
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load snpEff/4.3q

java -Xmx8G -jar $SNPEFF eff \
	-dataDir /isg/shared/databases/SnpEff/v4_3/data/ \
	hg19 na12878_q20_annot.vcf > na12878_q20_annot_snpEff.vcf


