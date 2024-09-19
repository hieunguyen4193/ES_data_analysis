gc()
rm(list = ls())

#####---------------------------------------------------------------------#####
##### libraries
#####---------------------------------------------------------------------#####
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

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

#####---------------------------------------------------------------------#####
##### input data: read the sub-cluster Myeloid_Basophils from integrated data all_sobj.integrated.rds
#####---------------------------------------------------------------------#####
outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
sub.cluster.id <- "Myeloid_Basophils"
s.obj.name <- "all_sobj.integrated.rds"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
path.to.08.output <- file.path(path.to.main.output, "08_output", sub.cluster.id)
dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(path.to.04.output, s.obj.name, "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", sub.cluster.id)))
s.obj <- subset(s.obj, name %in% c(sample1, sample2))

#####---------------------------------------------------------------------#####
##### remove TCR/BCR genes
#####---------------------------------------------------------------------#####
TR_genes_patterns <- c("Trav", "Traj", "Trac", "Trbv", "Trbd", "Trbj", "Trbc",
                       "Trgv", "Trgj", "Trgc", "Trdv", "Trdc", "Trdj") 
DefaultAssay(s.obj) <- "RNA"

all_genes <- row.names(s.obj)

genes.to.exclude <- unlist(lapply(all_genes, function(x){
  if (substr(x, 1, 4) %in% TR_genes_patterns){
    return(x)
  } else {
    return(NA)
  }
}))
genes.to.exclude <- subset(genes.to.exclude, is.na(genes.to.exclude) == FALSE)
genes.to.keep <- setdiff(row.names(s.obj), genes.to.exclude)

#####---------------------------------------------------------------------#####
##### add cluster annotation and mark cells that are in contamination ambient RNA clusters
#####---------------------------------------------------------------------#####
annotationdf <- read.csv(file.path(path.to.main.src, "Myeloid_Basophils_cluster_annotations.csv"))#
ambientRNA.clusters <- c(4, 5, 6, 8, 9, 10, 13, 17,  18, 19)
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(celltype = subset(annotationdf, annotationdf$Cluster == seurat_clusters)$CellType) %>%
  mutate(ambientRNA.cluster = ifelse(seurat_clusters %in% ambientRNA.clusters, "yes", "no")) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data),]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$celltype, col.name = "celltype")
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$ambientRNA.cluster, col.name = "ambientRNA.cluster")

#####---------------------------------------------------------------------#####
##### remove B cells and T cells contamination
#####---------------------------------------------------------------------#####
source(file.path(path.to.main.src, "run_integration.R"))
s.obj <- subset(s.obj, celltype %in% c("T cell contaminated", "B cell contaminated") == FALSE)

#####---------------------------------------------------------------------#####
##### extract cells that are  ambient RNA contaminated
#####---------------------------------------------------------------------#####
s.obj.ambientRNA <- subset(s.obj, ambientRNA.cluster == "yes")

#####---------------------------------------------------------------------#####
##### do re-integration for the new seurat objects
#####---------------------------------------------------------------------#####
if (file.exists(file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.rds", sub.cluster.id))) == FALSE){
  s.obj <- run_integration_after_subsetting(s.obj)
  saveRDS(s.obj, file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.rds", sub.cluster.id)))
  
} else {
  s.obj <- readRDS(file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.rds", sub.cluster.id)))
}

if (file.exists(file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.ambientRNAonly.rds", sub.cluster.id))) == FALSE){
  s.obj.ambientRNA <- run_integration_after_subsetting(s.obj.ambientRNA)
  saveRDS(s.obj.ambientRNA, file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.ambientRNAonly.rds", sub.cluster.id)))
  
} else {
  s.obj.ambientRNA <- readRDS(file.path(path.to.08.output, sprintf("%s.remove_contam_T_B_cells.ambientRNAonly.rds", sub.cluster.id)))
}