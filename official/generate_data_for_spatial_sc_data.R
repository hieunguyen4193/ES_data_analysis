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
path.to.save.cloupe.file <- file.path(path.to.main.output, "cloupe", integration.case, regression.mode, filter.mode)
dir.create(path.to.save.cloupe.file, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.10.output <- file.path(path.to.main.output, "10_output", integration.case, regression.mode, filter.mode)

path.to.s.obj <- file.path(path.to.10.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
s.obj <- readRDS(path.to.s.obj)

##### sub-sampling each clusters of the main seurat object
sampling.rate <- 0.75

cluster.barcodes <- list()
subsampling.cluster.barcodes <- list()

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")

set.seed(411)

all.subsampling.cells <- c()
for (cluster.id in unique(meta.data$cca.cluster.0.5)){
  cluster.barcodes[[cluster.id]] <- subset(meta.data, meta.data$cca.cluster.0.5 == cluster.id)$barcode
  subsampling.cluster.barcodes[[cluster.id]] <- sample(cluster.barcodes[[cluster.id]], round(sampling.rate * length(cluster.barcodes[[cluster.id]]) ))
  all.subsampling.cells <- c(all.subsampling.cells, 
                             subsampling.cluster.barcodes[[cluster.id]])
}

s.obj <- subset(s.obj, cells = all.subsampling.cells)
s.obj.no.reInt <- s.obj

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
cluster.resolution <- 0.5

DefaultAssay(s.obj) <- "RNA"
s.obj <- JoinLayers(s.obj)
s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                     save.RDS.s8 = TRUE,
                                                     path.to.output = path.to.save.cloupe.file,
                                                     use.sctransform = TRUE,
                                                     num.PCA = num.PCA,
                                                     num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                     num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                     cluster.resolution = cluster.resolution,
                                                     vars.to.regress = vars.to.regress)

##### CONVERT seurat object to cloupe file
if ("loupeR" %in% installed.packages() == FALSE){
  install.packages("hdf5r")
  install.packages("/media/hieunguyen/HD01/storage/offline_pkgs/loupeR_Linux.tar.gz", repos = NULL, type = "source")
}

library(loupeR)
loupeR::setup()
create_loupe_from_seurat(
  s.obj,
  output_dir = file.path(path.to.save.cloupe.file),
  output_name = sprintf("converted_cloupe_file.subSampling_%s.reIntegration.cloupe", sampling.rate),
  dedup_clusters = FALSE,
  executable_path = NULL,
  force = TRUE)

create_loupe_from_seurat(
  s.obj.no.reInt,
  output_dir = file.path(path.to.save.cloupe.file),
  output_name = sprintf("converted_cloupe_file.subSampling_%s.cloupe", sampling.rate),
  dedup_clusters = FALSE,
  executable_path = NULL,
  force = TRUE)

