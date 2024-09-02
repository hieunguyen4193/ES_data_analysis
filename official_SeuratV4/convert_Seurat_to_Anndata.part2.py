import scanpy as sc
import numpy as np
import os
import anndata2ri
import pathlib
import scvelo as scv
from scipy import io
import anndata#
import pandas as pd
from tqdm import tqdm
import argparse
import sys

# Activate the anndata2ri conversion between SingleCellExperiment and AnnData
anndata2ri.activate()

#Loading the rpy2 extension enables cell magic to be used
#This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython

sc.settings.verbosity = 3
# sc.logging.print_versions()

import warnings
warnings.filterwarnings("ignore")

outdir = "/home/hieunguyen/CRC1382/outdir"
all_datasets = ["full", "Myeloid_Basophils", "T_cells", "B_cells"]
PROJECT = "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
for dataset_name in tqdm(all_datasets):
    path_to_main_output = os.path.join(outdir, PROJECT, "data_analysis")
    path_to_seurat2anndata = os.path.join(path_to_main_output, "seurat2anndata")
    #####---------------------------------------------------------------------------------------------------------#####
    # load sparse matrix:
    sample_name = dataset_name
    X = io.mmread(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(sample_name)))
    
    # create anndata object
    adata = anndata.AnnData(X=X.transpose().tocsr())
    
    # load cell metadata:
    cell_meta = pd.read_csv(os.path.join(path_to_seurat2anndata, "metadata_{}.csv".format(sample_name)))
    
    # load gene names:
    with open(os.path.join(path_to_seurat2anndata, "gene_names_{}.csv".format(sample_name)), 'r') as f:
      gene_names = f.read().splitlines()
    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names
    
    # load dimensional reduction:
    pca = pd.read_csv(os.path.join(path_to_seurat2anndata, "pca_{}.csv".format(sample_name)))
    pca.index = adata.obs.index
    
    # set pca and umap
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
    
    # save dataset as anndata format
    print("writting file adata to {}".format(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(sample_name))))
    adata.write(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(sample_name)))
    
    # reload dataset
    adata = sc.read_h5ad(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(sample_name)))
