#####----------------------------------------------------------------------#####
##### Creating Single Cell References for Xenium Custom Panel Design from Seurat or AnnData
#####----------------------------------------------------------------------#####

gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))
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

##### update 07.01.2025
##### subset: keep only Hashtag1
s.obj <- subset(s.obj, HTO_classification == "Hashtag1-TotalSeqC")

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


##### annotation
annot <- list(
  `0` = "Naive T cells",
  `1` = "Fully differentiated effectors",
  `2` = "Mature B cells",
  `3` = "Immature neutrophils",
  `4` = "NK/ILC1s",
  `5` = "Pro B cells",
  `6` = "Monocytes/moMacs",
  `7` = "NKT cells",
  `8` = "Naive T cells 2",
  `9` = "Mature neutrophils",
  `10` = "Treg",
  `11` = "pDC",
  `12` = "Cycling T cells",
  `13` = "GMP",
  `14` = "Pre B cells",
  `15` = "DC",
  `16` = "Prdx1 high T cells",
  `17` = "Exhausted T cells",
  `18` = "Circualting neutrophils",
  `19` = "Kupffer cells",
  `20` = "Basophils",
  `21` = "Erythroids",
  `22` = "Migratory DC",
  `23` = "Other B cells")

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(celltype = annot[[cca.cluster.0.5]]) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$celltype, col.name = "celltype")
  
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

# DefaultAssay(s.obj) <- "RNA"
# s.obj <- JoinLayers(s.obj)
# 
# if (file.exists(file.path(path.to.save.cloupe.file, "s8_output", "EStange_20240411_reduced_RNAcontam_0.output.s8.rds")) == FALSE){
#   print("running integration ...")
#   s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
#                                                        save.RDS.s8 = TRUE,
#                                                        path.to.output = path.to.save.cloupe.file,
#                                                        use.sctransform = TRUE,
#                                                        num.PCA = num.PCA,
#                                                        num.PC.used.in.UMAP = num.PC.used.in.UMAP,
#                                                        num.PC.used.in.Clustering = num.PC.used.in.Clustering,
#                                                        cluster.resolution = cluster.resolution,
#                                                        vars.to.regress = vars.to.regress, 
#                                                        PROJECT = sprintf("%s_subsampling", PROJECT))
# } else {
#   print("File exists. reading in  ...")
#   s.obj.integrated <- readRDS(file.path(path.to.save.cloupe.file, "s8_output", "EStange_20240411_reduced_RNAcontam_0.output.s8.rds"))
# }

##### CONVERT 
source(
  "/home/hieunguyen/CRC1382/src_2023/EStange/official/convert_MEX_helper_functions.R"
)

if ("R.utils" %in% installed.packages() == FALSE){
  install.packages("R.utils")
}

# Define function
# savedir <- file.path(path.to.save.cloupe.file, "newVersion_02122024")
# seurat_obj <- s.obj.integrated
# annot.col.name <- "celltype"

run_save_mex <- function(seurat_obj, savedir, annot.col.name){
  # Run function
  print("running writeCounts ...")
  count.mat <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
  excluded.genes <- setdiff(row.names(count.mat), row.names(s.obj))
  print(dim(count.mat[setdiff(row.names(count.mat), excluded.genes), ]))
  writeCounts(
    savedir,
    count.mat[setdiff(row.names(count.mat), excluded.genes), ],
    gene.id = rownames(seurat_obj),
    gene.symbol = GetAssay(seurat_obj)@meta.features[["feature_name"]],
    feature.type = GetAssay(seurat_obj)@meta.features[["feature_type"]],
    barcodes = colnames(seurat_obj)
  )
  print("finished writeCounts")
  print(seurat_obj)
  print("start bundleOuputs ...")
  bundleOutputs(out_dir = savedir, data = seurat_obj, cell_type = annot.col.name)
  p <- DimPlot(object = seurat_obj, reduction = "cca_UMAP", label = TRUE, label.box = TRUE, group.by = annot.col.name)
  ggsave(plot = p, filename = "UMAP_annotated.svg", path = savedir, dpi = 300, width = 14, height = 10, device = "svg")
}

# run_save_mex(seurat_obj = s.obj.integrated, savedir = file.path(path.to.save.cloupe.file, "newVersion_025122024_reIntegrated"), annot.col.name = "celltype")
run_save_mex(seurat_obj = s.obj.no.reInt, savedir = file.path(path.to.save.cloupe.file, "newVersion_20250107"), annot.col.name = "celltype")
