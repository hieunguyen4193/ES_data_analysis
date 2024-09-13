gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"
chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25
num.PC.used.in.UMAP <- 25
my_random_seed <- 42

#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####
# integration.case <- params$integration.case
# regression.mode <- params$regression.mode
# filter.mode <- params$filter.mode

integration.case <- "remove_d4_LPS_SC5"
regression.mode <- "CC_differences"
filter.mode <- "nCount_and_BCR_TCRgenes"

PROJECT <- "EStange_20240411_reduced_RNAcontam_0"
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.10.output <- file.path(path.to.main.output, "10_output", integration.case, regression.mode, filter.mode)

path.to.s.obj <- file.path(path.to.10.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
s.obj <- readRDS(path.to.s.obj)

#### CONVERT seurat object to cloupe file
if ("loupeR" %in% installed.packages() == FALSE){
  install.packages("hdf5r")
  install.packages("/media/hieunguyen/HD01/storage/offline_pkgs/loupeR_Linux.tar.gz", repos = NULL, type = "source")
}

if (file.exists(file.path(path.to.07.output, sprintf("PROJECT_%s_%s_cloupe_converted_from_seurat", PROJECT, 
                                                     sub.cluster.idx))) == FALSE){
  library(loupeR)
  loupeR::setup()
  create_loupe_from_seurat(
    s.obj,
    output_dir = file.path(path.to.07.output),
    output_name = sprintf("PROJECT_%s_%s_cloupe_converted_from_seurat", PROJECT, 
                          sub.cluster.idx),
    dedup_clusters = FALSE,
    executable_path = NULL,
    force = TRUE)  
}

