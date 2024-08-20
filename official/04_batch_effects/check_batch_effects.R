#####----------------------------------------------------------------------#####
#
# trnguyen@ukaachen.de
#
#####----------------------------------------------------------------------#####

##### clean up #####
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
input.case <- "remove_d4_LPS"

PROJECT <- "EStange_20240411_SeuratV5"
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.03.output <- file.path(path.to.main.output, "03_output", input.case)
path.to.02.output <- file.path(path.to.main.output, "02_output")

path.to.04.output <- file.path(path.to.main.output, "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

s.obj.SC5 <- readRDS(file.path(path.to.02.output, "SC5", sprintf("GEX_sample_%s_seurat_object.rds", "SC5")))
s.obj.d10_SPF <- readRDS(file.path(path.to.02.output, "d10_SPF", sprintf("GEX_sample_%s_seurat_object.rds", "d10_SPF")))

s.obj.be <- merge(x = s.obj.SC5, s.obj.d10_SPF)
s.obj.be <- DietSeurat(s.obj.be)

my_random_seed <- 42
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
pca_reduction_name <- "RNA_PCA"
cluster.resolution <- 0.5
umap_reduction_name <- "RNA_UMAP"

s.obj.be <- SCTransform(s.obj.be, vars.to.regress = vars.to.regress, verbose = FALSE)
s.obj.be <- RunPCA(s.obj.be, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
s.obj.be <- RunUMAP(s.obj.be, reduction = pca_reduction_name, 
                 dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                 seed.use = my_random_seed, umap.method = "uwot")
# clustering 
s.obj.be <- FindNeighbors(s.obj.be, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
s.obj.be <- FindClusters(s.obj.be, resolution = cluster.resolution, random.seed = 0)

p <- DimPlot(object = s.obj.be, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name")
ggsave(plot = p, filename = "Sample_SC5_and_d10_SPF.svg", path = path.to.04.output, device = "svg", dpi = 300, width = 14, height = 10)
p <- DimPlot(object = s.obj.be, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, split.by = "name")
ggsave(plot = p, filename = "Sample_SC5_and_d10_SPF.split_by_name.svg", path = path.to.04.output, device = "svg", dpi = 300, width = 14, height = 10)
