#!/bin/bash
#SBATCH --job-name=index
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

module load bwa/0.7.17
bwa index -p chr20 chr20.fa
module unload bwa/0.7.17



module load samtools/1.7
samtools faidx chr20.fa


