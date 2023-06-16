#!/bin/bash 

#SBATCH --job-name=hhseq_merge
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err

chip_number_1="E100051605"
chip_number_2="E100057013"

lib="16"

str1=_

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/

mkdir $chip_number_1$str1$chip_number_2

cd $chip_number_1$str1$chip_number_2

mkdir $lib

cd $lib

mkdir rawdata

cd rawdata

zcat /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1/$lib/rawdata/H_R1.fq.gz /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_2/$lib/rawdata/H_R1.fq.gz > /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1$str1$chip_number_2/$lib/rawdata/H_R1.fq && gzip /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1$str1$chip_number_2/$lib/rawdata/H_R1.fq

zcat /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1/$lib/rawdata/H_R2.fq.gz /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_2/$lib/rawdata/H_R2.fq.gz > /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1$str1$chip_number_2/$lib/rawdata/H_R2.fq && gzip /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number_1$str1$chip_number_2/$lib/rawdata/H_R2.fq
