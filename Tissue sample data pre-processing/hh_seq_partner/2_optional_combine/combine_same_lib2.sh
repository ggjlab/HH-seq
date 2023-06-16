#!/bin/bash 

#SBATCH --job-name=hhseq_merge
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%J.out
#SBATCH --error=%J.err

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/E100041558_E100051718/79+80+13/processed/fastq_correction

zcat E100051718/H_R1_trimmed_modified_final.fastq.gz E100041558/H_R1_trimmed.fastq.gz > H_R1_trimmed.fastq && gzip H_R1_trimmed.fastq

zcat E100051718/H_R2_trimmed.fastq.gz E100041558/H_R2_trimmed.fastq.gz > H_R2_trimmed.fastq && gzip H_R2_trimmed.fastq