#!/bin/bash

script=cell_number_infer.sh

while read line
do

	chip_number=`echo $line | awk '{print $1}'`
	batch=`echo $line | awk '{print $2}'`
	index_type=`echo $line | awk '{print $3}'`
	species=`echo $line | awk '{print $4}'`
	cell_number_preset=`echo $line | awk '{print $5}'`
	
	sbatch $script $chip_number $batch $index_type $species $cell_number_preset

done < ./cell_number_preset.txt;
