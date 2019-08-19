#!/bin/bash
#SBATCH --job-name=sickle
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

module load sickle
sickle pe -t sanger \
	-f ../01_raw_data/NA12878_20_paired_1.fq -r ../01_raw_data/NA12878_20_paired_2.fq \
	-o trimmed_NA12878_20_paired_1.fq -p trimmed_NA12878_20_paired_2.fq \
	-l 45 \
	-q 25 \
	-s singles_NA12878_20_paired_2.fq  

module load fastqc
fastqc -o ../03_fastqc/ trimmed_NA12878_20_paired_1.fq trimmed_NA12878_20_paired_2.fq 


