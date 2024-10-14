gc()
rm(list = ls())

library(monocle)
library(Seurat)
library(stringr)
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.15.output <- file.path(path.to.main.output, "15_output", "input_monocle2_objs")
dir.create(path.to.15.output, showWarnings = FALSE, recursive = TRUE)

all.s.obj <- read.csv("/home/hieunguyen/CRC1382/src_2023/EStange/official/15_monocle2/list_of_objects_to_run_monocle2.csv", sep = ",")

for (path.to.s.obj in all.s.obj$path){
  print(sprintf("Working on %s", path.to.s.obj))
  split.name <- str_split(path.to.s.obj, "/")[[1]]
  output.index <- split.name[[9]]
  integration.case <- split.name[[10]]
  regresison.mode <- split.name[[11]]
  filter.mode <- split.name[[12]]
  sub.cluster <- split.name[[13]]
  
  if (file.exists(file.path(path.to.15.output, 
                            sprintf("from_%s", output.index),
                            integration.case,
                            regresison.mode,
                            filter.mode,
                            sub.cluster,
                            sprintf("%s.rds", sub.cluster))) == FALSE){
    s.obj <- readRDS(path.to.s.obj)
    
    data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
    
    pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
    
    fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fd)
    
    library(monocle)
    monocle.obj <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
    dir.create(file.path(path.to.15.output, 
                         sprintf("from_%s", output.index),
                         integration.case,
                         regresison.mode,
                         filter.mode,
                         sub.cluster),
               showWarnings = FALSE, recursive = TRUE)
    saveRDS(monocle.obj, file.path(path.to.15.output, 
                                   sprintf("from_%s", output.index),
                                   integration.case,
                                   regresison.mode,
                                   filter.mode,
                                   sub.cluster,
                                   sprintf("%s.rds", sub.cluster)))  
  }
}

