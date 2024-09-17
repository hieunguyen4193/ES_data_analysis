gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.samplesheet <- file.path(path.to.main.src, "sampleSheets_for_DGE_and_CellChat")
dir.create(path.to.save.samplesheet, showWarnings = FALSE, recursive = TRUE)

##### get s.obj objects from the 03 output
all.03.results <- Sys.glob(file.path(path.to.main.output, "03_output", "*", "s8_output", "*.rds"))
tmp.03.output <- data.frame(
  output_index = c("03_output"),
  integration.case = to_vec( for(item in all.03.results) str_split(item, "/")[[1]][[10]]),
  path = all.03.results
)

tmp.03.output[["num.clusters"]] <- unlist(
  lapply(tmp.03.output$path, function(x){
    tmp.s.obj <- readRDS(x)
    return(length(unique(tmp.s.obj$cca.cluster.0.5)))
  })
)
writexl::write_xlsx(tmp.03.output, file.path(path.to.save.samplesheet, "SampleSheet_03_output.xlsx"))

##### get s.obj objects from the 10 output
all.10.results <- Sys.glob(file.path(path.to.main.output, "10_output", "*", "*", "*", "s8_output", "*.rds"))
tmp.10.output <- data.frame(
  output_index = c("10_output"),
  integration.case = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[10]]),
  regression.mode = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[11]]),
  filter.mode = to_vec( for(item in all.10.results) str_split(item, "/")[[1]][[12]]),
  path = all.10.results
)
tmp.10.output[["num.clusters"]] <- unlist(
  lapply(tmp.10.output$path, function(x){
    tmp.s.obj <- readRDS(x)
    return(length(unique(tmp.s.obj$cca.cluster.0.5)))
  })
)
writexl::write_xlsx(tmp.10.output, file.path(path.to.save.samplesheet, "SampleSheet_10_output.xlsx"))
