#!/bin/bash 

#SBATCH --job-name=hhseq_n_cell
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err

#Contact: Jikai Shao & Xiunan Fang

########################################################################################################################
##############################	Presets
########################################################################################################################

chip_number=$1 # Chip number
batch=$2 # Sample index
index_type=$3 # Index type ('preindex' or 'no_preindex')
species=$4 # Species ('human', 'mouse' or 'hm_mixed')
cell_number_preset=$5

##########	Human & Mouse mixed genome/gtf

hm_mixed_genome=/public/home/guogjgroup/ggj/jikai/hm_mix_from_yongcheng/Mix_Human_Mus_Chu/STAR_index/
hm_mixed_gtf=/public/home/guogjgroup/ggj/jikai/hm_mix_from_yongcheng/Mix_Human_Mus_Chu/GRCh38_GRCm39.105.gtf

##########	Human genome/gtf

human_genome=/public/home/guogjgroup/ggj/jikai/STAR_2.7.10a_hg38/GenomeDir/
human_gtf=/public/home/guogjgroup/ggj/jikai/STAR_2.7.10a_hg38/Homo_sapiens.GRCh38.gtf

##########	Mouse genome/gtf

mouse_genome=/public/home/guogjgroup/ggj/jikai/STAR_2.7.10a_mm10/GenomeDir/
mouse_gtf=/public/home/guogjgroup/ggj/tools/STAR_Reference_Mouse/Mus_musculus.GRCm38.88.gtf

########################################################################################################################
##############################	Pre-processing
########################################################################################################################

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$batch/

mkdir processed

cd rawdata

cp H_R1.fq.gz H_R2.fq.gz ../processed/

cd ../processed/

mkdir fastq_correction

mv H_R1.fq.gz H_R2.fq.gz fastq_correction/

cd fastq_correction

# Sort R1 & R2 based on their 5' adapter
#bbduk2.sh in=H_R1.fq.gz outm=H_R1_linker2.fastq fliteral=GGAGTTGGAGTGAGTGGATGAGTGATG k=27 hdist=3 restrictleft=27 -Xmx80g

#bbduk2.sh in=H_R2.fq.gz outm=H_R2_linker2.fastq fliteral=GGAGTTGGAGTGAGTGGATGAGTGATG k=27 hdist=3 restrictleft=27 -Xmx80g

#bbduk2.sh in=H_R1.fq.gz outm=H_R1_linker4.fastq fliteral=GTGAGTGATGGTTGAGGATGTGTGGAGATA k=30 hdist=3 restrictleft=30 -Xmx80g

#bbduk2.sh in=H_R2.fq.gz outm=H_R2_linker4.fastq fliteral=GTGAGTGATGGTTGAGGATGTGTGGAGATA k=30 hdist=3 restrictleft=30 -Xmx80g

# Generate new R1 & R2
#cat H_R1_linker4.fastq H_R2_linker4.fastq > R1.fastq
#cat H_R1_linker2.fastq H_R2_linker2.fastq > R2.fastq

#/public/home/guogjgroup/ggj/tools/bbmap/repair.sh in=R1.fastq in2=R2.fastq out=R1_repaired.fastq out2=R2_repaired.fastq

# Sort R1 & R2 based on their 5' adapter
#bbduk2.sh in=R1_repaired.fastq in2=R2_repaired.fastq outm=R1_repaired_linker2.fastq outm2=R2_repaired_linker2.fastq  fliteral=GGAGTTGGAGTGAGTGGATGAGTGATG skipr1=t k=27 hdist=3 restrictleft=27 -Xmx80g

#bbduk2.sh in=R1_repaired.fastq in2=R2_repaired.fastq outm=R1_repaired_linker4.fastq outm2=R2_repaired_linker4.fastq  fliteral=GTGAGTGATGGTTGAGGATGTGTGGAGATA skipr1=t k=30 hdist=3 restrictleft=30 -Xmx80g

# Generate final R1 & R2
#cat R1_repaired_linker2.fastq R2_repaired_linker4.fastq > R1_final.fastq
#cat R2_repaired_linker2.fastq R1_repaired_linker4.fastq > R2_final.fastq

# Trim 5' adapters: R1 30bp & R2 27bp
# SOAPnuke: trim some bp of the read's head and tail, they means: (PE type:read1's head and tail and read2's head and tail [0,0,0,0]; SE type:read head and tail [0,0])
#echo -e "trim=30,0,27,0" >./soapnuke.config
#SOAPnuke filter -1 R1_final.fastq -2 R2_final.fastq -o ./ -C H_R1P.fq.gz -D H_R2P.fq.gz -J -c ./soapnuke.config

# Trim polyA
cutadapt -a AAAAAAAAAAAAAAA --minimum-length=30 --pair-filter=any -o H_R2_trimmed.fastq.gz -p H_R1_trimmed.fastq.gz H_R2.fq.gz H_R1.fq.gz > trimming_report.txt

rm H_R1.fq.gz H_R2.fq.gz

cd ..

########################################################################################################################
##############################	Whitelist-Set cell number
########################################################################################################################

mkdir barcode_extraction

cd barcode_extraction

# Automatic Identification of potential cell barcodes
# R1 Barcode(5'):|CB1(10bp)|TGGT(mismatches <= 2)|CB2(10bp)|GAGA(mismatches <= 2)|CB3(10bp)|UMI(8bp)
# R2 Barcode(preindex)(5'):CB4(10bp)
# Setting cell number cutoff to extract barcodes of real cell

if [ "$index_type"x = "preindex"x ]

then

	# Barcode = CB1+CB2+CB3+CB4
	umi_tools whitelist --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --bc-pattern2="(?P<cell_2>.{10}).*" --stdout whitelist_trimmed_extracted.txt --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz -L whitelist.log --ignore-read-pair-suffixes --plot-prefix=whitelist --filtered-out reads_filtered1.fastq --filtered-out2 reads_filtered2.fastq  --set-cell-number $cell_number_preset

else

	# Barcode = CB1+CB2+CB3
	umi_tools whitelist --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --stdout whitelist_trimmed_extracted.txt --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz -L whitelist.log --ignore-read-pair-suffixes --plot-prefix=whitelist --filtered-out reads_filtered1.fastq --filtered-out2 reads_filtered2.fastq  --set-cell-number $cell_number_preset


fi

cp whitelist_cell_barcode_counts.png /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/other/cell_number_preset/$batch.png

cd ..


