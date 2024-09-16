#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV5"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.13.output <- file.path(outdir, PROJECT, "data_analysis", "13_output")

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path.to.save.samplesheet <- file.path(path.to.main.src, "sampleSheets_for_DGE_and_CellChat")

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "13_DGE"
path.to.rmd <- file.path(path.to.main.src, src.dir, "13_DGE_analysis.Rmd")
output_dir <- file.path(path.to.save.html, "13_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
 `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output.xlsx")),
 `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output.xlsx"))
)

input.samplesheet <- samplesheets[["03_output"]]

path.to.save.html <- file.path(output_dir, "03_output")

i = 1
integration.case <- input.samplesheet[i, ][["integration.case"]]
path.to.s.obj <- input.samplesheet[i, ][["path"]]
output.file.name <- sprintf("%s_%s_%s_vs_%s.html", str_replace(basename(path.to.rmd), ".Rmd"), integration.case, sample1, sample2)
path.to.save.DGE.output <- file.path(path.to.13.output, integration.case, sprintf("%s_%s", sample1, sample2))
dir.create(path.to.save.DGE.output, showWarnings = FALSE, recursive = TRUE)

sample1 <- "adult_SPF"
sample2 <- "d7_SPF"

rmarkdown::render(path.to.rmd, 
                  params = list(
                    sample1 = sample1, 
                    sample2 = sample2, 
                    path.to.s.obj = path.to.s.obj,
                    path.to.save.output = path.to.save.DGE.output
                  ),
                  output_file = output.file.name,
                  output_dir = path.to.save.html)
