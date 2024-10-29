gc()
rm(list = ls())

source("/home/hieunguyen/CRC1382/src_2023/EStange/official/15_monocle2/15_monocle2_helper_functions.R")
my_random_seed <- 42

library(monocle)
library(Seurat)
library(stringr)
if ("svglite" %in% installed.packages() == FALSE){
  local({r <- getOption("repos")
  r["CRAN"] <- "http://cran.r-project.org"
  options(repos=r)})
  install.packages("svglite")
}
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.monocle.objs <- file.path(path.to.main.output, "15_output", "input_monocle2_objs")

all.monocle.objs <- Sys.glob(file.path(path.to.monocle.objs, "*/*/*/*/*/*.rds"))

for (path.to.monocle in all.monocle.objs[[3]]){
  print(sprintf("working on %s", path.to.monocle))
  split.name <- str_split(path.to.monocle, "/")[[1]]
  output.index <- split.name[[10]]
  integration.case <- split.name[[12]]
  regresison.mode <- split.name[[13]]
  filter.mode <- split.name[[14]]
  sub.cluster <- split.name[[15]] %>% str_replace(".rds", "")
  
  path.to.15.output <- file.path(path.to.main.output, 
                                 "15_output", 
                                 "monocle2_output",
                                 output.index,
                                 integration.case,
                                 regresison.mode,
                                 filter.mode,
                                 sub.cluster)
  dir.create(path.to.15.output, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.15.output, "finished.csv")) == FALSE){
    if (file.exists(file.path(path.to.15.output, "monocle_obj.rds")) == FALSE){
      print(sprintf("running monocle2 on dataset %s", sub.cluster))
      monocle.obj <- readRDS(path.to.monocle)
      monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.15.output)
    } else {
      print("monocle result exists, reading in...")
      monocle.obj <- readRDS(file.path(path.to.15.output, "monocle_obj.rds"))
    }
    
    ##### plot cell trajectory, color by seurat clusters
    p <- plot_cell_trajectory(monocle.obj, color_by = "seurat_clusters")
    ggsave(plot = p, filename = sprintf("cell_trajectory_%s.seurat_clsuters.svg", sub.cluster), path = path.to.15.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    ##### plot cell trajectory, color by monocle2 states
    p <- plot_cell_trajectory(monocle.obj, color_by = "State")
    ggsave(plot = p, filename = sprintf("cell_trajectory_%s.State.svg", sub.cluster), path = path.to.15.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    ##### plot cell trajectory, color by pseudotime
    p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
    ggsave(plot = p, filename = sprintf("cell_trajectory_%s.pseudotime.svg", sub.cluster), path = path.to.15.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    ##### plot cell trajectory, color by pseudotime
    monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
    p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
    ggsave(plot = p, filename = sprintf("cell_trajectory_%s.rev_Pseudotime.svg", sub.cluster), path = path.to.15.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    ##### save monocle data to csv file
    monocledf <- data.frame(
      barcode = colnames(monocle.obj),
      state = monocle.obj$State,
      pseudotime = monocle.obj$Pseudotime
    )
    monocle.reversedf <- data.frame(
      barcode = colnames(monocle.obj.reverse),
      state = monocle.obj.reverse$State,
      pseudotime = monocle.obj.reverse$Pseudotime
    )
    write.csv(monocledf, file.path(path.to.15.output, "monocledf.csv"))
    write.csv(monocle.reversedf, file.path(path.to.15.output, "monocledf.rev.csv"))
    
    write.csv(data.frame(status = c("finished")), file.path(path.to.15.output, "finished.csv"))
  } else {
    print(sprintf("Finished %s", path.to.15.output))
  }
}
