import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from pandas.testing import assert_frame_equal
from pysankey2 import Sankey
from scipy import stats
from pysankey2.datasets import load_fruits
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from scipy.spatial.distance import pdist
import random as rd

cancer_type = 'Glioma'
patient_list = ['Glioma_p1', 'Glioma_p2', 'Glioma_p3', 'Glioma_p4', 'Glioma_p5']

cnv_score_neighbor_overall = pd.DataFrame()

for patient_id in patient_list:
    
    # Automatically assign cancer_type

    if ('HCC' in patient_id) or ('ICC' in patient_id):

        cancer_type = 'HCC_ICC'

    elif ('ESCA' in patient_id):

        cancer_type = 'ESCA'

    elif ('LUAD' in patient_id):

        cancer_type = 'LUAD'

    elif ('COAD' in patient_id) or ('READ' in patient_id):

        cancer_type = 'COAD_READ'

    elif ('BRCA' in patient_id):

        cancer_type = 'BRCA'

    elif ('STAD' in patient_id):

        cancer_type = 'STAD'

    elif ('Glioma' in patient_id):

        cancer_type = 'Glioma'

    elif ('AML' in patient_id):

        cancer_type = 'AML'

    #gene_loc = pd.read_csv('/mnt/hdd1/jikai/cnv_test/gene_loc.tsv', sep = '\t', header=None)

    # Loading Clone

    # Loading inferCNV clone
    infer_cnv_clone = pd.read_csv('/mnt/hdd1/jikai/hhseq_processing_latest/inferCNV_per_cancer_type/' + cancer_type 
                                  + '_all/inferCNV/' + patient_id + 
                                  '_default_test/pyramidinal_101/infercnv_reshape.observation_groupings.txt', 
                                  sep = ' ')

    # Loading Numbat Clone
    #cnv_clone = pd.read_csv('/mnt/hdd1/jikai/numbat/' + patient_id + '/' + patient_id + '_cnv_clone_sorted.csv', 
    #                        index_col = 0)

    # Loading CNV Score

    respath = "/mnt/hdd1/jikai/hhseq_processing_latest/inferCNV_per_cancer_type/" + cancer_type + "_all/inferCNV/" + patient_id + "_default_test/pyramidinal_101/"

    # Loading CNV Score: Raw
    cnv_score_ob = pd.read_csv(respath + "infercnv_reshape.observations.txt",sep=' ')
    cnv_score_ob = cnv_score_ob.T.copy()

    # Loading CNV Score: bin
    #cnv_score_ob_cut = pd.read_csv(respath + 'cnv_score_ob_cut.csv', index_col=0)

    # Loading adata
    adata = sc.read('/mnt/hdd1/jikai/hhseq_processing_latest/annotation_per_patient/' + patient_id + '/' + 
                    patient_id + '_annotation.h5ad')

    '''
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html

    connectivitiessparse matrix of dtype float32.
    Weighted adjacency matrix of the neighborhood graph of data points. Weights should be interpreted as 
    connectivities.

    distancessparse matrix of dtype float32.
    Instead of decaying weights, this stores distances for each pair of neighbors.
    '''

    neighbor_df = adata.obsp['connectivities']
    neighbor_df = pd.DataFrame(neighbor_df.todense())
    neighbor_df.index = adata.obs.index.tolist()
    neighbor_df.columns = adata.obs.index.tolist()
    
    # Additional Settings
    n_pseudo_cell = 9
    cnv_score_used = cnv_score_ob.copy()
    cnv_clone_used = infer_cnv_clone.copy()
    clone_column = 'Dendrogram Group'
    
    cnv_score_neighbor_all = pd.DataFrame(index=cnv_score_used.columns.tolist())
    neighbor_df = neighbor_df.loc[cnv_score_used.index.tolist(), cnv_score_used.index.tolist()]

    for each_cell in tqdm(cnv_score_used.index.tolist()):

        each_cell_neighbor_subset = neighbor_df[[each_cell]]
        each_cell_neighbor_subset = each_cell_neighbor_subset[each_cell_neighbor_subset[each_cell] != 0]
        each_cell_neighbor_subset = pd.DataFrame(each_cell_neighbor_subset[each_cell].sort_values(ascending = False))
        each_cell_neighbor_subset = each_cell_neighbor_subset[:n_pseudo_cell-1]

        neighbor_cell_list = each_cell_neighbor_subset.index.tolist()
        neighbor_cell_list.append(each_cell)
        cnv_score_subset = cnv_score_used.loc[neighbor_cell_list, :]
        cnv_score_subset_neighbor = pd.DataFrame(cnv_score_subset.mean())
        cnv_score_subset_neighbor.columns = [each_cell]
        cnv_score_subset_neighbor = cnv_score_subset_neighbor.loc[cnv_score_neighbor_all.index.tolist(), :]
        cnv_score_neighbor_all[each_cell] = cnv_score_subset_neighbor[each_cell]

    cnv_score_neighbor_all = cnv_score_neighbor_all.T
    
    cnv_score_neighbor_overall = pd.concat([cnv_score_neighbor_overall, cnv_score_neighbor_all])

cnv_score_neighbor_overall = cnv_score_neighbor_overall.round(4)
cnv_score_neighbor_overall.to_csv(cancer_type + '_CNV_neighbor_all.csv')

