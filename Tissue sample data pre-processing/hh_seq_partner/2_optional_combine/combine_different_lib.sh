#!/bin/bash 

#SBATCH --job-name=hhseq_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%J.out
#SBATCH --error=%J.err

chip_number="E120220616"

str1="86"
str2="87"

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/

mkdir $str1+$str2

cd $str1+$str2

mkdir rawdata

cd rawdata

zcat /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1/rawdata/H_R1.fq.gz /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str2/rawdata/H_R1.fq.gz > /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1+$str2/rawdata/H_R1.fq && gzip /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1+$str2/rawdata/H_R1.fq

zcat /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1/rawdata/H_R2.fq.gz /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str2/rawdata/H_R2.fq.gz > /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1+$str2/rawdata/H_R2.fq && gzip /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$str1+$str2/rawdata/H_R2.fq