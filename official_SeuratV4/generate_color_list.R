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

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

all.datasets <- c("full", "Myeloid_Basophils", "T_cells", "B_cells")

for (sub.cluster.id in all.datasets){
  print(sprintf("working on dataset %s", sub.cluster.id))
  path.to.save.colors <- file.path(outdir, PROJECT, "colors", sub.cluster.id)
  dir.create(path.to.save.colors, showWarnings = FALSE, recursive = TRUE)
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
  num.clusters <- length(unique(s.obj$seurat_clusters))
  library(scales)
  colors <- hue_pal()(num.clusters)
  write.csv(data.frame(color = colors), file.path(path.to.save.colors, "colordf.csv"))
}