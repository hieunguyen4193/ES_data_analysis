import scanpy as sc
import numpy as np
import os
import anndata2ri
import pathlib
from scipy import io
import anndata#
import pandas as pd
from tqdm import tqdm
import argparse
import sys

# Activate the anndata2ri conversion between SingleCellExperiment and AnnData
# anndata2ri.activate()

#Loading the rpy2 extension enables cell magic to be used
#This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython

sc.settings.verbosity = 3
# sc.logging.print_versions()

import warnings
warnings.filterwarnings("ignore")

outdir = "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT = "EStange_20240411_reduced_RNAcontam_0"

path_to_main_src = "/home/hieunguyen/CRC1382/src_2023/EStange/official"
samplesheet = pd.read_excel(os.path.join(path_to_main_src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))

for row_i in range(len(samplesheet)):
    output_index = samplesheet.iloc[row_i]['output_index']
    input_info = []
    input_info_cols = ["integration.case", "regression.mode", "filter.mode", "sub.cluster.id"]

    for c in input_info_cols:
        if not pd.isna(samplesheet.iloc[row_i][c]):
            input_info.append(samplesheet.iloc[row_i][c])

    path_to_input_s_obj = samplesheet.iloc[row_i]['path']
    path_to_main_output = os.path.join(outdir, PROJECT)

    path_to_18_output = os.path.join(path_to_main_output, "18_output", f"from_{output_index}")
    path_to_seurat2anndata = os.path.join(path_to_18_output, "seurat2anndata", "/".join(input_info))
        
    print(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(PROJECT)))
    X = io.mmread(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(PROJECT)))

    # create anndata object
    adata = anndata.AnnData(X=X.transpose().tocsr())

    # load cell metadata:
    cell_meta = pd.read_csv(os.path.join(path_to_seurat2anndata, "metadata_{}.csv".format(PROJECT)))

    # load gene names:
    with open(os.path.join(path_to_seurat2anndata, "gene_names_{}.csv".format(PROJECT)), 'r') as f:
        gene_names = f.read().splitlines()

    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names

    # load dimensional reduction:
    pca = pd.read_csv(os.path.join(path_to_seurat2anndata, "pca_{}.csv".format(PROJECT)))
    pca.index = adata.obs.index

    # set pca and umap
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

    # save dataset as anndata format
    adata.write(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(PROJECT)))
