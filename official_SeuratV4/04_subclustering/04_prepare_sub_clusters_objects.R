##### clean up #####
# gc()
# rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
PROJECT <- "230605_220928_Stange"
run <- "230605_220928_Stange"

force_rerun <- FALSE

analysis.round <- "1st"

chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25

path.to.src <- "/home/hieunguyen/CRC1382/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(path.to.src, "s8_integration_and_clustering.R"))

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir/EStange"

sub.clusters <- hash()

sub.clusters[["all_sobj_remove_d4_LPS.integrated.rds"]] <- list(Lymphoid = c(0, 1, 4, 5, 10, 13, 18, 16),
                                                                B_cells = c(2, 8, 15, 17, 22, 25),
                                                                Myeloid_Neutrophil = c(3, 6, 7, 9, 11, 12, 14, 19, 24, 20, 21, 23, 26, 27),
                                                                Myeloid_only =  c(11, 12, 14, 19, 24, 20, 21))

sub.clusters[["all_sobj.integrated.rds"]] <- list(Myeloid_MHCII = c(8, 10, 12, 15, 22, 27))

# loop through all integrated R.objects
for (s.obj.name in c("all_sobj_remove_d4_LPS.integrated.rds", "all_sobj.integrated.rds")){
  sub.cluster.idxs <- sub.clusters[[s.obj.name]]
  path.to.main.output <- file.path("/media/hieunguyen/CRC1382H/CRC1382/outdir/EStange/CRC1382_EStange_projects", run, "data_analysis")
  
  path.to.02.output <- file.path(path.to.main.output, "02_output")
  
  path.to.03.output <- file.path(path.to.main.output, "03_output")
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.04.output <- file.path(path.to.main.output, "04_output", s.obj.name)
  dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
  
  dir.create(file.path(path.to.04.output, "sub_cluster_objects"), showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.04.output, "FINISHED_STEP04.csv")) == FALSE){
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
      saveRDS(object = sub.cluster.obj, file = file.path(path.to.04.output, "sub_cluster_objects", sprintf("sub_cluster_%s.rds", group)))
    }
    saveRDS(sub.clusterdf, file.path(path.to.04.output, "sub_cluster_objects", "convert_cluster_to_celltype_table.rds"))
    write.csv(data.frame(data = c("Finished generating subset of cells at step 04")), file.path(path.to.04.output, "FINISHED_STEP04.csv"))  
  }

}

write.csv(data.frame(data = c("Finished generating subset of cells at step 04")), file.path(path.to.04.output, "FINISHED_STEP04_ALL_CASES.csv"))
