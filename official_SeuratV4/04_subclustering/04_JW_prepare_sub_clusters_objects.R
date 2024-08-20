# gc()
# rm(list = ls())

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

path.to.src <- "/home/hieunguyen/CRC1382/src/src_pipeline/scRNA_GEX_pipeline/processes_src"

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

sub.clusters <- hash()

sub.clusters[["all_sobj.integrated.rds"]] <- list(Myeloid_Basophils = c(5, 3, 13, 18, 23, 22, 7, 10, 8, 24, 15, 27, 12),
                                                  B_cells = c(21, 16, 9, 19, 26, 2, 25),
                                                  T_cells = c(17, 20, 6, 14, 4, 11, 1, 0),
                                                  subset_231031 = c(8, 10, 7, 22, 13, 3, 5))

# loop through all integrated R.objects
for (s.obj.name in c("all_sobj.integrated.rds")){
  sub.cluster.idxs <- sub.clusters[[s.obj.name]]
  path.to.04.output <- file.path(path.to.main.output, "04_output", s.obj.name)
  dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.04.output, "sub_cluster_objects"), showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.04.output, "FINISHED_STEP04_JW.csv")) == FALSE){
    s.obj <- readRDS(file.path(path.to.03.output, s.obj.name))
    
    sub.clusterdf <- data.frame()
    for (group in names(sub.cluster.idxs)){
      tmpdf <- data.frame(sub.cluster.idxs[[group]])
      tmpdf$celltype <- group
      colnames(tmpdf) <- c("cluster", "celltype")
      sub.clusterdf <- rbind(sub.clusterdf, tmpdf)
    }
    
    for (group in names(sub.cluster.idxs)){
      sub.cluster.obj <- subset(s.obj, seurat_clusters %in% sub.cluster.idxs[[group]])
      saveRDS(object = sub.cluster.obj, file = file.path(path.to.04.output, "sub_cluster_objects", sprintf("sub_cluster_%s_JW.rds", group)))
    }
    saveRDS(sub.clusterdf, file.path(path.to.04.output, "sub_cluster_objects", "convert_cluster_to_celltype_table.rds"))
    write.csv(data.frame(data = c("Finished generating subset of cells at step 04")), file.path(path.to.04.output, "FINISHED_STEP04_JW.csv"))  
  }

}

write.csv(data.frame(data = c("Finished generating subset of cells at step 04")), file.path(path.to.04.output, "FINISHED_STEP04_ALL_CASES_JW.csv"))
