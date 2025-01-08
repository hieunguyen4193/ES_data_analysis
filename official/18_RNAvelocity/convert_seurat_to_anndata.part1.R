##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))

for (row_i in seq(1, nrow(samplesheet))){
  output.index <- samplesheet[row_i, ]$output_index
  input.info <- c()
  input.info.cols <- c("integration.case", 
                       "regression.mode", 
                       "filter.mode", 
                       "sub.cluster.id")
  
  for (c in input.info.cols){
    if (is.na(samplesheet[row_i, ][[c]]) == FALSE){
      input.info <- c(input.info, samplesheet[row_i, ][[c]])
    }
  }
  path.to.input.s.obj <- samplesheet[row_i, ]$path
  path.to.main.output <- file.path(outdir, PROJECT)
  
  path.to.18.output <- file.path(path.to.main.output, "18_output", sprintf("from_%s", output.index))
  path.to.seurat2anndata <- file.path(path.to.18.output, "seurat2anndata", paste0(input.info, collapse = "/"))
  
  dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', PROJECT))) == FALSE){
    print("generating ...")
    
    s.obj <- readRDS(path.to.input.s.obj)
    
    s.obj$barcode <- colnames(s.obj)
    
    s.obj$UMAP_1 <- s.obj@reductions$cca_UMAP@cell.embeddings[,1]
    s.obj$UMAP_2 <- s.obj@reductions$cca_UMAP@cell.embeddings[,2]
    
    write.csv(s.obj@reductions$cca_UMAP@cell.embeddings, 
              file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', PROJECT)), 
              quote=F, 
              row.names=F)
    
    write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', PROJECT)), quote=F, row.names=F)
    
    # write expression counts matrix
    library(Matrix)
    counts_matrix <- GetAssayData(s.obj, assay='SCT', slot='data')
    writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', PROJECT)))
    
    # write gene names
    write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', PROJECT)),
                 quote=F,row.names=F,col.names=F)
    
    coldf <- data.frame(cluster = unique(s.obj$cca.cluster.0.5),
                        color = hue_pal()(length(unique(s.obj$cca.cluster.0.5))))
    write.csv(coldf, file.path(path.to.seurat2anndata, sprintf('colordf_%s.csv', PROJECT)))
  } else {
    print(sprintf("Seurat2Anndata for project %s already exists at %s", PROJECT, path.to.seurat2anndata))
  }
}