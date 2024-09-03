gc()
rm(list = ls())

# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

my_random_seed <- 42
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
library(monocle)

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
source(file.path(path.to.project.src, "helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"

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

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
s.obj.name <- "all_sobj.integrated.rds"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")

path.to.monocle2.input <- file.path(path.to.main.output, "monocle2_inputs")
dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)


for (sub.cluster.id in c("full", "Myeloid_Basophils", "T_cells", "B_cells")){
  monocle.obj <- readRDS(file.path(path.to.monocle2.input, sprintf("%s.rds", sub.cluster.id)))  
  path.to.06.output <- file.path(path.to.main.output, "06_output", sub.cluster.id)
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.06.output, "monocle_obj.rds")) == FALSE){
    print(sprintf("running monocle2 on dataset %s", sub.cluster.id))
    monocle.obj <- readRDS(file.path(path.to.monocle2.input, sprintf("%s.rds", sub.cluster.id)))
    monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.06.output)
  } else {
    print("monocle result exists, reading in...")
    monocle.obj <- readRDS(file.path(path.to.06.output, "monocle_obj.rds"))
  }
  
  ##### plot cell trajectory, color by seurat clusters
  p <- plot_cell_trajectory(monocle.obj, color_by = "seurat_clusters")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.seurat_clsuters.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by monocle2 states
  p <- plot_cell_trajectory(monocle.obj, color_by = "State")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.State.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.pseudotime.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
  
  ##### plot cell trajectory, color by pseudotime
  monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
  p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
  ggsave(plot = p, filename = sprintf("cell_trajectory_%s.rev_Pseudotime.svg", sub.cluster.id), path = path.to.06.output, device = "svg", dpi = 300, width = 14, height = 10)
  
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
  write.csv(monocledf, file.path(path.to.06.output, "monocledf.csv"))
  write.csv(monocle.reversedf, file.path(path.to.06.output, "monocledf.rev.csv"))
}
