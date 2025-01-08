##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
library(monocle)

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))
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
  
  dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.monocle2.input, 
                            sprintf("%s.monocle2.rds", PROJECT))) == FALSE){
    print(sprintf("File %s does not exists, generating ...",
                  file.path(path.to.monocle2.input, 
                            sprintf("%s.monocle2.rds", PROJECT))))
    s.obj <- readRDS(path.to.input.s.obj)
    
    data <- GetAssayData(s.obj, slot = "data", assay = "SCT")
    
    pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
    
    fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fd)
    
    monocle.obj <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
    saveRDS(monocle.obj, file.path(path.to.monocle2.input, 
                                   sprintf("%s.monocle2.rds", PROJECT)))
  } else {
    print(sprintf("File %s exists", file.path(path.to.monocle2.input, 
                                              sprintf("%s.monocle2.rds", PROJECT))))
  }
}