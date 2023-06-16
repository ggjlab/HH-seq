#! /usr/bin/Rscript

library('Seurat')
library('infercnv')
reticulate::use_condaenv("/mnt/hdd1/tools/anaconda3/", required = TRUE)
options(scipen = 100)

window_used <- 101
smooth_used <- "pyramidinal"

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = 'raw.tsv',
                                    annotations_file = 'anno.tsv',
                                    gene_order_file = '/mnt/hdd1/jikai/cnv_test/gene_loc.tsv',
                                    ref_group_names = c('Endothelial Cell', 'Fibroblast', 'Macrophage')
                                    )

infercnv_obj <- suppressWarnings(infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             cluster_by_groups = FALSE,
                             analysis_mode = "subclusters",
                             tumor_subcluster_partition_method = "leiden",
                             k_obs_groups = 7, 
                             smooth_method = smooth_used,
                             window_length = window_used,
                             num_threads = 20,
                             max_centered_threshold = "auto",
                             denoise = TRUE,
                             HMM = FALSE,
                             out_dir = paste(smooth_used, "_", window_used, "/", sep = "")))

plot_cnv(infercnv_obj,
                 out_dir=paste(smooth_used, "_", window_used, "/", sep = ""),
                 title="inferCNV",
                 output_filename="infercnv_reshape",
                 output_format="png",
                 write_expr_matrix=FALSE,
                 cluster_by_groups = FALSE,
                 k_obs_groups = 7,
                 png_res=300,
                 useRaster=FALSE,
                 dynamic_resize = 0.02)

