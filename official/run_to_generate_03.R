#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "03_data_integration"
path.to.rmd <- file.path(path.to.main.src, src.dir, "03_analysis_on_integrated_dataset.Rmd")
output_dir <- file.path(path.to.save.html, "03_output")

for (input.case in c("remove_d4_LPS",
                     "all.samples",
                     "remove_d4_LPS_SC5")){
  output_file <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", input.case))
  
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(input.case = input.case))}  
}
