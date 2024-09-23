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
path.to.10.output <- file.path(path.to.main.output, "10_output", sub.cluster.id)
dir.create(path.to.10.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(path.to.08.output, "Myeloid_Basophils.remove_contam_T_B_cells.rds"))
s.obj <- subset(s.obj, name %in% c(sample1, sample2))

all.s.obj <- list(
  d4_LPS = subset(s.obj, name == "d4_LPS"),
  d4_SPF = subset(s.obj, name == "d4_SPF")
)

cluster.markers <- list()
if (file.exists(file.path(path.to.10.output, "sample_cluster_markers.rds")) == FALSE){
  for (sample.id in names(all.s.obj)){
    tmp <- FindAllMarkers(object = all.s.obj[[sample.id]], assay = "RNA", test.use = "wilcox")
    cluster.markers[[sample.id]] <- subset(tmp, tmp$p_val_adj <= 0.05 & tmp$avg_log2FC > 0)    
  }
  saveRDS(cluster.markers, file.path(path.to.10.output, "sample_cluster_markers.rds"))
} else {
  cluster.markers <- readRDS(file.path(path.to.10.output, "sample_cluster_markers.rds"))
}

available.clusters <- intersect(
  unique(cluster.markers$d4_LPS$cluster),
  unique(cluster.markers$d4_SPF$cluster)
)

for (cluster.id in available.clusters){
  dir.create(file.path(path.to.10.output, sprintf("compare_clusterDEG_cluster%s", cluster.id)), showWarnings = FALSE, recursive = TRUE)
  de.gene.LPS <- subset(cluster.markers$d4_LPS, cluster.markers$d4_LPS$cluster == cluster.id)$gene
  de.gene.SPF <- subset(cluster.markers$d4_SPF, cluster.markers$d4_SPF$cluster == cluster.id)$gene
  
  intersect.genes <- intersect(de.gene.LPS, de.gene.SPF)
  in.LPS.only <- setdiff(de.gene.LPS, de.gene.SPF)
  in.SPF.only <- setdiff(de.gene.SPF, de.gene.LPS)
  writexl::write_xlsx(data.frame(intersect.genes = intersect.genes),
                      file.path(path.to.10.output, sprintf("compare_clusterDEG_cluster%s", cluster.id), "DEG_in_both_LPS_SPF.xlsx"))
  writexl::write_xlsx(data.frame(in.LPS.only = in.LPS.only),
                      file.path(path.to.10.output, sprintf("compare_clusterDEG_cluster%s", cluster.id), "DEG_in_LPS_only.xlsx"))
  writexl::write_xlsx(data.frame(in.SPF.only = in.SPF.only),
                      file.path(path.to.10.output, sprintf("compare_clusterDEG_cluster%s", cluster.id), "DEG_in_SPF_only.xlsx"))
}

