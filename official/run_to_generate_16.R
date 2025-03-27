#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs", "16_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "16_monocle3"
path.to.rmd <- file.path(path.to.main.src, src.dir, "16_monocle3_analysis.Rmd")

path.to.16.output <- file.path(path.to.main.output, "16_output")
dir.create(path.to.16.output, showWarnings = FALSE, recursive = TRUE)

samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))
samplelist <- read.csv(file.path(path.to.main.src, "13_DGE", "sample_comparision_list.csv"))

library(argparse)
parser <- ArgumentParser()

parser$add_argument("-i", "--row_i", action="store",
                    help="Full name of the input project/dataset name")

args <- parser$parse_args()

row_i <- args$row_i

# for (row_i in seq(1, nrow(samplesheet))){
    path.to.s.obj <- samplesheet[row_i, ]$path
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
    
    path.to.save.output <- file.path(
      path.to.16.output, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/")
    )
    
    output_file <- sprintf("16_monocle3_analysis.html")
    
    output_dir <- file.path(
      path.to.save.html, 
      sprintf("from_%s", samplesheet[row_i, ]$output_index),
      paste0(input.info, collapse = "/")
    )
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    excluded.clusters <- NULL
    use_partition <- FALSE
    reduction.name <- "cca_UMAP"
    cluster.name <- "cca.cluster.0.5"
    
    input.params <- list(
      path.to.s.obj = path.to.s.obj,
      path.to.save.output = path.to.save.output,
      use_partition = use_partition,
      excluded.clusters = excluded.clusters,
      reduction.name = reduction.name, 
      cluster.name = cluster.name
    )
    for (n in names(input.params)){
      print(sprintf("Input param %s = %s", n, input.params[[n]]))
    }
    
    if (file.exists(file.path(output_dir, output_file)) == FALSE){
      print(sprintf("Generate new html file %s...", file.path(output_dir, output_file) ))
      rmarkdown::render(input = path.to.rmd,
                        output_file = output_file,
                        output_dir = output_dir,
                        params = input.params)
    } else {
      print(sprintf("html file %s exists...", file.path(output_dir, output_file) ))
    }
  # }
