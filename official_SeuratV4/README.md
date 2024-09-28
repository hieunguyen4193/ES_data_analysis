# Data analysis: Two samples `d4_LPS` vs `d4_SPF`

## Downstream analysis pipeline

From the output of `CellRanger`, we run the pipeline (backbone `Seurat`) to perform some Quality controls, downstream analysis and data integration. Run `run_pipeline_RNAcontam_thresholds.R`. 

## Data analysis pipeline

### 01: Hashtag antibodies quality control and cell hashtag assignment

Run `run_to_generate_01.R`. This script calls all *.Rmd* files in the folder `01_QC_hashtag_antibodies_data`. 

### 02: Perform downstream analysis for each sample

Further downstream analysis was performed when calling `run_to_generate_02.R`. This scripts call `*.Rmd` files in `02_preliminary_analysis`. Note that, we removed **doublet** hashtag cells and keep **negative** hashtag cells. Differential gene expression analysis between clusters were also conducted, which might help identifying cell populations. 

### 03: Perform data integration for all samples (not only `d4_LPS` and `d4_SPF`)
Run `run_to_generate_03.R`. 

### 04: Sub-clustering based on samples and based on cell types

- Run `04_subclustering.Rmd`: in this analysis, we take the integrated data object generated in `03` and extract selected clusters (See the image `sub_clustering.png` attached). Output objects are saved at `04_output/all_sobj.integrated.rds/preprocessed_subcluster_obj/preprocessed_subclusters_<sub_cluster_name>.rds`.

```r
sub.clusters <- hash()
sub.clusters[["all_sobj.integrated.rds"]] <- list(
  Myeloid_Basophils = c(20, 5, 3, 8, 17, 16, 25, 26, 9, 7, 12, 28, 11),
  B_cells = c(2, 21, 27, 29, 15, 10, 22),
  T_cells = c(14, 1, 6, 13, 4, 24, 0))
```

- Run `04_subset_d4_samples.Rmd`: this also takes the integrated data object, but extract cells from samples `d4_LPS` and `d4_SPF` only. Output object is saved at `04_output/all_sobj.integrated_d4_LPS_d4_SPF_only.rds`.

### 05: Perform differential gene expression analysis between `d4_LPS` vs `d4_SPF`
Run `run_to_generate_05_*.R` (two versions: for 04 output and for 08 output) to perform differential gene expression analysis. This function calls the Rmarkdown files `05_DGE_analysis.pseudoBulk_cluster.Rmd` and `05_DGE_analysis.Rmd`. 

### Further downstream analysis based on Ambient RNA level of the cells

Run `prepare_data_for_further_analysis_with_ambientRNA_cells.R` to prepare the data: select cells based on Ambient RNA contamination levels. 


