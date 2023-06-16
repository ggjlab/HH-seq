#!/bin/bash 

#SBATCH --job-name=hh_seq_processing
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=16
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
##############################	Whitelist-Set cell number
########################################################################################################################

cd /public/home/guogjgroup/ggj/rawdata/MATQ_Pan_Cancer/human_pan_cancer/$chip_number/$batch/

cd processed

#mkdir barcode_extraction

cd barcode_extraction

# Automatic Identification of potential cell barcodes
# R1 Barcode(5'):|CB1(10bp)|TGGT(mismatches <= 2)|CB2(10bp)|GAGA(mismatches <= 2)|CB3(10bp)|UMI(8bp)
# R2 Barcode(preindex)(5'):CB4(10bp)
# Setting cell number cutoff to extract barcodes of real cell

if [ "$index_type"x = "preindex"x ]

then

	# Barcode = CB1+CB2+CB3+CB4
	umi_tools whitelist --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --bc-pattern2="(?P<cell_2>.{10}).*" --stdout whitelist_trimmed_extracted.txt --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz -L whitelist.log --ignore-read-pair-suffixes --plot-prefix=whitelist --filtered-out reads_filtered1.fastq --filtered-out2 reads_filtered2.fastq  --set-cell-number $cell_number_preset

	# Extract barcodes and UMIS and add to read names
	umi_tools extract --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --bc-pattern2="(?P<cell_2>.{10}).*" --stdout H_R1_trimmed_extracted.fastq --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz --read2-out=H_R2_trimmed_extracted.fastq --whitelist=whitelist_trimmed_extracted.txt -L extract.log --ignore-read-pair-suffixes

else

	# Barcode = CB1+CB2+CB3
	umi_tools whitelist --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --stdout whitelist_trimmed_extracted.txt --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz -L whitelist.log --ignore-read-pair-suffixes --plot-prefix=whitelist --filtered-out reads_filtered1.fastq --filtered-out2 reads_filtered2.fastq  --set-cell-number $cell_number_preset

	# Extract barcodes and UMIS and add to read names
	umi_tools extract --stdin ../fastq_correction/H_R1_trimmed.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{30})(?P<umi_1>.{8}).*" --stdout H_R1_trimmed_extracted.fastq --read2-in ../fastq_correction/H_R2_trimmed.fastq.gz --read2-out=H_R2_trimmed_extracted.fastq --whitelist=whitelist_trimmed_extracted.txt -L extract.log --ignore-read-pair-suffixes

fi

cd ..

########################################################################################################################
##############################	Mapping & DGE generation
########################################################################################################################

mkdir star_mapping

cd star_mapping

if [ "$species"x = "hm_mixed"x ]

then 

	genome_file=$hm_mixed_genome
	gtf_file=$hm_mixed_gtf

elif [ "$species"x = "human"x ]

then

	genome_file=$human_genome
	gtf_file=$human_gtf

else

	genome_file=$mouse_genome
	gtf_file=$mouse_gtf

fi

# STAR Mapping (2 pass & generation of unmapped reads)
STAR --runThreadN 6 --runMode alignReads --genomeDir $genome_file --readFilesIn ../barcode_extraction/H_R2_trimmed_extracted.fastq --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 41143265264 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --chimOutJunctionFormat 1 --chimSegmentMin 5 --outReadsUnmapped Fastx --twopassMode Basic --limitOutSJcollapsed 2000000

# extract uniquely mapped reads
# The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10*log10(1- 1/Nmap)) for multi-mapping reads.
samtools view -b -q 250 Aligned.sortedByCoord.out.bam | samtools sort - > star.bam

cd ..

mkdir feature_count

cd feature_count

# Assign gene feature, with strandness
featureCounts -t gene -g gene_name -a $gtf_file -o gene_assigned.bw25113 -R BAM ../star_mapping/star.bam -T 10 -s 2 --extraAttributes gene_biotype
samtools sort star.bam.featureCounts.bam -o assigned.sorted.bam
samtools index assigned.sorted.bam

cd ..

mkdir dge

cd dge 

# dge generation
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I ../feature_count/assigned.sorted.bam -S counts.tsv

# Summary barcode reads
path=`pwd`
input_folder=$path
output_folder=$path
fastqfile=../barcode_extraction/H_R2_trimmed_extracted.fastq
python /public/home/guogjgroup/ggj/tools/preindex_MATQ/Barcode_reads.py $input_folder $output_folder $fastqfile > barcode_read_count.tsv
sort -k2 -nr barcode_read_count.tsv > barcode_read_count_sort.tsv

python /public/home/guogjgroup/ggj/tools/preindex_MATQ/Python_script/UMI_Gene_count_xin.py counts.tsv barcode_read_count_sort.tsv UMI_Gene_count.tsv
sed -i '1d' UMI_Gene_count.tsv
sort -k3 -nr UMI_Gene_count.tsv > UMI_Gene_count_sorted.tsv

cd ..

########################################################################################################################
##############################	Exon_only Generation
########################################################################################################################

cd feature_count

mkdir exon_only

cd exon_only

# Assign gene feature, with strandness
featureCounts -t exon -g gene_name -a $gtf_file -o exon_assigned.bw25113 -R BAM ../../star_mapping/star.bam -T 10 -s 2 -J --extraAttributes gene_biotype
samtools sort star.bam.featureCounts.bam -o assigned.sorted.bam
samtools index assigned.sorted.bam

cd ../../dge

mkdir exon_only

cd exon_only

# dge generation
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I ../../feature_count/exon_only/assigned.sorted.bam -S counts.tsv

# Summary barcode reads
python /public/home/guogjgroup/ggj/tools/preindex_MATQ/Python_script/UMI_Gene_count_xin.py counts.tsv ../barcode_read_count_sort.tsv UMI_Gene_count.tsv
sed -i '1d' UMI_Gene_count.tsv
sort -k3 -nr UMI_Gene_count.tsv > UMI_Gene_count_sorted.tsv

cd ../..

########################################################################################################################
##############################	Human mouse mixture processing
########################################################################################################################

if [ "$species"x = "hm_mixed"x ]

then 

	mkdir hm_mixed

	cp ../feature_count/assigned.sorted.bam hm_mixed/

	cd hm_mixed

	/public/home/guogjgroup/ggj/tools/Drop-seq_tools-2.5.1/FilterBam I=./assigned.sorted.bam O=./mouse.bam REF_SOFT_MATCHED_RETAINED=GRCm39
	samtools index mouse.bam

	/public/home/guogjgroup/ggj/tools/Drop-seq_tools-2.5.1/FilterBam I=./assigned.sorted.bam O=./human.bam REF_SOFT_MATCHED_RETAINED=hg38
	samtools index human.bam

	mkdir mouse
	mv mouse.bam mouse.bam.bai mouse/

	mkdir human
	mv human.bam human.bam.bai human/

	cp ../barcode_read_count_sort.tsv human/
	cp ../barcode_read_count_sort.tsv mouse/

	# Mouse
	cd mouse/

	umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I mouse.bam -S counts.tsv
	python /public/home/guogjgroup/ggj/tools/preindex_MATQ/Python_script/UMI_Gene_count_xin.py counts.tsv barcode_read_count_sort.tsv UMI_Gene_count.tsv

	sed -i '1d' UMI_Gene_count.tsv
	sort -k3 -nr UMI_Gene_count.tsv > m_UMI_Gene_count_sorted.tsv

	# Huamn
	cd ../human/

	umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I human.bam -S counts.tsv
	python /public/home/guogjgroup/ggj/tools/preindex_MATQ/Python_script/UMI_Gene_count_xin.py counts.tsv barcode_read_count_sort.tsv UMI_Gene_count.tsv
	
	sed -i '1d' UMI_Gene_count.tsv
	sort -k3 -nr UMI_Gene_count.tsv > h_UMI_Gene_count_sorted.tsv

fi



