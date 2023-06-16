#!/bin/bash 

#SBATCH --job-name=hhseq_split
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%J.out
#SBATCH --error=%J.err

# Presets
chip_number=E100066004
batch=(41
	42
	43
	44
	45
	46
	47
	48
	62
	64
	65
	)

str1=_L01_read_1.fq.gz
str2=_L01_read_2.fq.gz

str3=_1.fq.gz
str4=_2.fq.gz
str5=_L01_

# Splitting Barcode
cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/L01/

/public/home/guogjgroup/ggj/tools/splitBarcode-master/V0.1.6_release/linux/splitBarcode \
/public/home/guogjgroup/ggj/tools/splitBarcode-master/index_HUADA.txt $chip_number$str1 -2 $chip_number$str2 -o OUTPUT -b 202 10 2 -r

# Sorting Outputs
for i in ${batch[@]}

do

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/

mkdir $i

cd $i

mkdir rawdata

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/L01/OUTPUT/

mv $chip_number$str5$i$str3 /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$i/rawdata/
mv $chip_number$str5$i$str4 /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$i/rawdata/

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$i/rawdata/

mv $chip_number$str5$i$str3 H_R1.fq.gz
mv $chip_number$str5$i$str4 H_R2.fq.gz

done

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/

mv L01 other

cd other

mkdir cell_number_preset


