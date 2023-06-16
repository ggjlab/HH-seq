#!/bin/sh
#SBATCH --job-name=SCENIC_1
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=%J.out
#SBATCH --error=%J.err

patient_id=SCENIC_Epithelium_subset_1
n_core=32

str1=_adj.tsv
str2=_raw_counts.csv

pyscenic grn --num_workers $n_core --output /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str1 --method grnboost2 /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str2 /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/cisTarget_resources/tf_lists/allTFs_hg38.txt

str1=_adj.tsv
str2=_raw_counts.csv
str3=_reg.csv

pyscenic ctx /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str1 /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/cisTarget_resources/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/cisTarget_resources/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/cisTarget_resources/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str2 --mode "dask_multiprocessing" --output /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str3 --num_workers $n_core

str1=_adj.tsv
str2=_raw_counts.csv
str3=_reg.csv
str4=_SCENIC.csv

str5=_all_cell_raw_counts.csv

pyscenic aucell /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str5 /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str3 --output /public/home/guogjgroup/ggj/matq_analysis/pan_cancer/SCENIC/$patient_id/$patient_id$str4 --num_workers $n_core