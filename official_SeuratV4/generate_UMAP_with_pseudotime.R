##### clean up #####
gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
s.obj.name <- "all_sobj.integrated.rds"

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

all.datasets <- c("full", "Myeloid_Basophils", "T_cells", "B_cells")

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

for (sub.cluster.id in all.datasets){
  print(sprintf("working on dataset %s", sub.cluster.id))
  path.to.save.colors <- file.path(outdir, PROJECT, "colors", sub.cluster.id)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")  
  path.to.03.output <- file.path(path.to.main.output, "03_output")
  path.to.04.output <- file.path(path.to.main.output, "04_output")
  path.to.05.output <- file.path(path.to.main.output, "05_output", sub.cluster.id)
  path.to.06.output <- file.path(path.to.main.output, "06_output", sub.cluster.id)
  
  if (sub.cluster.id == "full"){
    s.obj <- readRDS(file.path(path.to.03.output, s.obj.name))
  } else {
    s.obj <- readRDS(file.path(path.to.04.output, s.obj.name, "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", sub.cluster.id)))  
  }
  s.obj <- subset(s.obj, name %in% c(sample1, sample2))
  
  monocledf <- read.csv(file.path(path.to.06.output, "monocledf.csv"))
  monocledf.rev <- read.csv(file.path(path.to.06.output, "monocledf.rev.csv"))
  colnames(monocledf.rev) <- c("X", "barcode", "state", "rev.pseudotime")
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") 
  meta.data <- merge(meta.data, subset(monocledf, select = c(barcode, pseudotime)), by.x = "barcode", by.y = "barcode")
  meta.data <- merge(meta.data, subset(monocledf.rev, select = c(barcode, rev.pseudotime)), by.x = "barcode", by.y = "barcode")
  
  meta.data <- meta.data %>% column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$pseudotime, col.name = "pseudotime") 
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$rev.pseudotime, col.name = "rev.pseudotime") 
  
  umap.reduction.name <- "INTE_UMAP"
  umap.pseudotime <- FeaturePlot(object = s.obj, reduction = umap.reduction.name, label = TRUE, features = c("pseudotime"))
  umap.rev.pseudotime <- FeaturePlot(object = s.obj, reduction = umap.reduction.name, label = TRUE, features = c("rev.pseudotime"))
  
  ggsave(plot = umap.pseudotime, filename = sprintf("%s_umap_pseudotime.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
  ggsave(plot = umap.rev.pseudotime, filename = sprintf("%s_umap_rev_pseudotime.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
}

