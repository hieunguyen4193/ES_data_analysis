gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"

chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

if ("heatmaply" %in% installed.packages() == FALSE){
  install.packages("heatmaply")
}
path.to.src <- "/home/hieunguyen/CRC1382/src/src_pipeline/scRNA_GEX_pipeline/processes_src"

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
sub.cluster.id <- "Myeloid_Basophils"
sub.cluster.id <- "T_cells"
sub.cluster.id <- "B_cells"
s.obj.name <- "all_sobj.integrated.rds"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
path.to.05.output <- file.path(path.to.main.output, "05_output", sub.cluster.id)
path.to.seurat2anndata <- file.path(path.to.main.output, "seurat2anndata")
dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)

for (sub.cluster.id in c("full", "Myeloid_Basophils", "T_cells", "B_cells")){
  if (file.exists(file.path(path.to.seurat2anndata, sprintf('finished_%s.csv', sub.cluster.id))) == FALSE){
    print(sprintf("working on the case %s", sub.cluster.id))
    if (sub.cluster.id == "full"){
      s.obj <- readRDS(file.path(path.to.03.output, s.obj.name))
    } else {
      s.obj <- readRDS(file.path(path.to.04.output, s.obj.name, "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", sub.cluster.id)))  
    }
    s.obj <- subset(s.obj, name %in% c(sample1, sample2))
    s.obj$barcode <- colnames(s.obj)
    
    s.obj$UMAP_1 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,1]
    s.obj$UMAP_2 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,2]
    write.csv(s.obj@reductions$INTE_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', sub.cluster.id)), quote=F, row.names=F)
    write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', sub.cluster.id)), quote=F, row.names=F)
    
    # write expression counts matrix
    library(Matrix)
    counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='data')
    writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', sub.cluster.id)))
    
    # write gene names
    write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', sub.cluster.id)),
                 quote=F,row.names=F,col.names=F)
    write.csv(data.frame(status = c("finished converting for the dataset")), file.path(path.to.seurat2anndata, sprintf('finished_%s.csv', sub.cluster.id)))    
  }
}