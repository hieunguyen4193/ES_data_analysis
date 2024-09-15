gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

##### get s.obj objects from the 03 output
all.03.results <- Sys.glob(file.path(path.to.main.output, "03_output", "*", "s8_output", "*.rds"))
tmp.03.output <- data.frame(
  output_index = c("03_output"),
  integration.case = to_vec( for(item in all.03.results) str_split(item, "/")[[1]][[9]]),
  path = all.03.results
)

##### get s.obj objects from the 10 output
all.10.results <- Sys.glob(file.path(path.to.main.output, "10_output", "*", "*", "*", "s8_output", "*.rds"))
tmp.10.output <- data.frame(
  output_index = c("10_output"),
  integration.case = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[10]]),
  regression.mode = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[11]]),
  filter.mode = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[12]]),
  path = all.10.results
)
