#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs", "14_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "14_CellChat"
path.to.rmd <- file.path(path.to.main.src, src.dir, "14_CellChat_diff_analysis_2_samples.Rmd")

path.to.14.output <- file.path(path.to.main.output, "14_output")
dir.create(path.to.14.output, showWarnings = FALSE, recursive = TRUE)

samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))
samplelist <- read.csv(file.path(path.to.main.src, "13_DGE", "sample_comparision_list.csv"))


for (row_i in seq(1, nrow(samplesheet))){
  for (row_j in seq(1, nrow(samplelist))){
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
    
    sample1 <- samplelist[row_j, ]$sample1
    sample2 <- samplelist[row_j, ]$sample2
    path.to.s.obj <- samplesheet[row_i, ]$path
    filter10cells <- "Filter10"
    
    path.to.cellchat1 <- file.path(
      path.to.14.output, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/"), 
      sample1,
      sprintf("CellChat_object.%s.%s.rds", sample1, filter10cells))
    
    path.to.cellchat2 <- file.path(
      path.to.14.output, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/"), 
      sample2,
      sprintf("CellChat_object.%s.%s.rds", sample2, filter10cells))
    
    path.to.save.output <- file.path(
      path.to.14.output, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/"),
      sprintf("%s_vs_%s", sample1, sample2)
    )
    
    output_file <- sprintf("CellChat_%s_vs_%s.html", sample1, sample2)
    
    output_dir <- file.path(
      path.to.save.html, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/"),
      sprintf("%s_vs_%s", sample1, sample2)
    )
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (file.exists(file.path(output_dir, output_file)) == FALSE){
      print(sprintf("Generate new html file %s...", file.path(output_dir, output_file) ))
      rmarkdown::render(input = path.to.rmd,
                        output_file = output_file,
                        output_dir = output_dir,
                        params = list(
                          path.to.s.obj = path.to.s.obj,
                          sample1 = sample1,
                          path.to.cellchat1 = path.to.cellchat1,
                          sample2 = sample2,
                          path.to.cellchat2 = path.to.cellchat2,
                          path.to.save.output = path.to.save.output,
                          filter10cells = filter10cells))
    } else {
      print(sprintf("html file %s exists...", file.path(output_dir, output_file) ))
    }
    
  }
}