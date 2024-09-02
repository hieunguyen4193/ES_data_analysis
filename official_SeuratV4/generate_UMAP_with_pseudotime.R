##### clean up #####
gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/migratory_ccr7_data_analysis"

source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

all.datasets <- c("220907_FH",
                  "GSM5764259", 
                  "integrate_GSE192742_LIVER",
                  "230228_FH",
                  "GSM5764288",            
                  "GSM5764245",
                  "gutcellatlas_myeloid",
                  "220907_FH_cDC1",
                  "220907_FH_cDC2"
)

outdir <- "/media/hieunguyen/HNSD01/outdir"

config.version <- "v0.1"
output.version <- "20240828"
PROJECT <- "FHager_datasets"

input.config <- config.params[[config.version]]

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
for (orig.dataset in all.datasets){
  print(sprintf("working on dataset %s", orig.dataset))
  dataset.name <- sprintf("%s_%s", orig.dataset, config.version)
  path.to.save.colors <- file.path(outdir, PROJECT, output.version, "colors", dataset.name)
  dir.create(path.to.save.colors, showWarnings = FALSE, recursive = TRUE)
  path.to.main.input <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis", "01_output", sprintf("%s.rds", dataset.name))
  path.to.main.output <- file.path(outdir, PROJECT, output.version, dataset.name, "data_analysis")
  path.to.02.output <- file.path(path.to.main.output, "02_output_monocle_2.12.0_v1")
  
  s.obj <- readRDS(path.to.main.input)
  monocledf <- read.csv(file.path(path.to.02.output, "monocledf.csv"))
  monocledf.rev <- read.csv(file.path(path.to.02.output, "monocledf.rev.csv"))
  colnames(monocledf.rev) <- c("X", "barcode", "state", "rev.pseudotime")
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") 
  meta.data <- merge(meta.data, subset(monocledf, select = c(barcode, pseudotime)), by.x = "barcode", by.y = "barcode")
  meta.data <- merge(meta.data, subset(monocledf.rev, select = c(barcode, rev.pseudotime)), by.x = "barcode", by.y = "barcode")
  
  meta.data <- meta.data %>% column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$pseudotime, col.name = "pseudotime") 
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$rev.pseudotime, col.name = "rev.pseudotime") 
  
  if (orig.dataset == "integrate_GSE192742_LIVER"){
    umap.reduction.name <- "INTE_UMAP"
  } else {
    umap.reduction.name <- "RNA_UMAP"
  }
  umap.pseudotime <- FeaturePlot(object = s.obj, reduction = umap.reduction.name, label = TRUE, features = c("pseudotime"))
  umap.rev.pseudotime <- FeaturePlot(object = s.obj, reduction = umap.reduction.name, label = TRUE, features = c("rev.pseudotime"))
  
  ggsave(plot = umap.pseudotime, filename = sprintf("%s_umap_pseudotime.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
  ggsave(plot = umap.rev.pseudotime, filename = sprintf("%s_umap_rev_pseudotime.svg", dataset.name), path = path.to.02.output, device = "svg", dpi = 300, width = 14, height = 10)
}

