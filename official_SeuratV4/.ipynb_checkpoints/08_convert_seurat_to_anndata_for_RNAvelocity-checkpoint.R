gc()
rm(list = ls())
scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
source(file.path(path.to.project.src, "monocle2_helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"

my_random_seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4"
s.obj.name <- "all_sobj.integrated.rds"
sub.cluster.id <- "subset_231031"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
path.to.08.output <- file.path(path.to.main.output, "08_output")
dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)

path.to.seurat2anndata <- path.to.08.output
path.to.sobj <- file.path(path.to.04.output, s.obj.name, "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", sub.cluster.id))
s.obj <- readRDS(path.to.sobj)
object.name <- str_replace(basename(path.to.sobj), ".rds", "")
s.obj <- subset(s.obj, name %in% c(sample1, sample2))
s.obj$barcode <- colnames(s.obj)
  
s.obj$UMAP_1 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,1]
s.obj$UMAP_2 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,2]

# write dimesnionality reduction matrix, in this example s.obj.case pca matrix
write.csv(s.obj@reductions$INTE_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', object.name)), quote=F, row.names=F)
write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', object.name)), quote=F, row.names=F)
  
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', object.name)))
  
# write gene names
write.table(data.frame('gene'=rownames(counts_matrix)), 
            file = file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', object.name)),
            quote = F, 
            row.names = F, 
            col.names = F)


