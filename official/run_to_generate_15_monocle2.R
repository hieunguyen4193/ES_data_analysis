##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)
library(dplyr)

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))
source(file.path(path.to.main.src, "15_monocle2", "monocle2_helper_functions.R"))

for (row_i in seq(1, nrow(samplesheet))){
  output.index <- samplesheet[row_i, ]$output_index
  input.info <- c()
  input.info.cols <- c("integration.case", 
                       "regression.mode", 
                       "filter.mode", 
                       "sub.cluster.id")
  
  for (c in input.info.cols){
    if (is.na(samplesheet[row_i, ][[c]]) == FALSE){
      input.info <- c(input.info, samplesheet[row_i, ][[c]])
    }
  }
  
  path.to.input.s.obj <- samplesheet[row_i, ]$path
  path.to.main.output <- file.path(outdir, PROJECT)
  
  path.to.15.output <- file.path(path.to.main.output, "15_output", sprintf("from_%s", output.index))
  path.to.monocle2.input <- file.path(path.to.15.output, "monocle2_input", paste0(input.info, collapse = "/"))
  path.to.monocle2.output <- file.path(path.to.15.output, "monocle2_output", paste0(input.info, collapse = "/"))
  
  monocle.obj <- readRDS(file.path(path.to.monocle2.input, sprintf("%s.monocle2.rds", PROJECT)))
  
  if (file.exists(file.path(path.to.monocle2.output, "monocle_obj.rds")) == FALSE){
    print(sprintf("running monocle2 on dataset %s", PROJECT))
    monocle.obj <- readRDS(path.to.monocle.obj)
    monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.monocle2.output)
  } else {
    print("monocle result exists, reading in...")
    monocle.obj <- readRDS(file.path(path.to.monocle2.output, "monocle_obj.rds"))
  }
  
  ##### plot cell trajectory, color by seurat clusters
  p <- plot_cell_trajectory(monocle.obj, color_by = "cca.cluster.0.5")
  ggsave(plot = p, 
         filename = sprintf("cell_trajectory_%s.seurat_clsuters.svg", PROJECT), 
         path = path.to.monocle2.output, 
         device = "svg", 
         dpi = 300, 
         width = 14, 
         height = 10)
  
  ##### plot cell trajectory, color by monocle2 states
  p <- plot_cell_trajectory(monocle.obj, color_by = "State")
  ggsave(plot = p, 
         filename = sprintf("cell_trajectory_%s.State.svg", PROJECT), 
         path = path.to.monocle2.output, 
         device = "svg", 
         dpi = 300, 
         width = 14, 
         height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.pseudotime.svg", PROJECT), 
         path = path.to.monocle2.output, 
         device = "svg", 
         dpi = 300, 
         width = 14, 
         height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
  p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
  ggsave(plot = p, 
         filename = sprintf("cell_trajectory_%s.rev_Pseudotime.svg", PROJECT), 
         path = path.to.monocle2.output, 
         device = "svg", 
         dpi = 300, 
         width = 14, 
         height = 10)
  
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
  write.csv(monocledf, file.path(path.to.monocle2.output, "monocledf.csv"))
  write.csv(monocle.reversedf, file.path(path.to.monocle2.output, "monocledf.rev.csv"))
}