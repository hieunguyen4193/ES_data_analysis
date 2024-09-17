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

path.to.13.output <- file.path(outdir, PROJECT, "data_analysis", "13_output")

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path.to.save.samplesheet <- file.path(path.to.main.src, "sampleSheets_for_DGE_and_CellChat")

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "14_CellChat"
path.to.rmd <- file.path(path.to.main.src, src.dir, "14_CellChat_general_analysis.Rmd")
output_dir <- file.path(path.to.save.html, "13_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
  `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output.xlsx")),
  `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output.xlsx"))
)

output.index <- "03_output"
input.sheet <- samplesheets[[output.index]]
i <- 1
sample.id <- "d7_GF"
path.to.s.obj <- input.sheet[i, ][["path"]]
filter10cells <- "Filter10"


