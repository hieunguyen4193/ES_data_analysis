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

src.dir <- "13_DGE"
path.to.rmd <- file.path(path.to.main.src, src.dir, "13_DGE_analysis.Rmd")
output.dir <- file.path(path.to.save.html, "13_output")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
 `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output.xlsx")),
 `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output.xlsx"))
)

comparison.samplesheet <- read.csv(file.path(path.to.main.src, src.dir, "sample_comparision_list.csv"))

for (output.index in names(samplesheets)){
  input.samplesheet <- samplesheets[[output.index]]
  for (row_i in seq(1, nrow(input.samplesheet))){
    for (j in seq(1, nrow(comparison.samplesheet))){
      sample1 <- comparison.samplesheet[row_i, ][["sample1"]]
      sample2 <- comparison.samplesheet[row_i, ][["sample2"]]
      
      if (output.index == "03_output"){
        integration.case <- input.samplesheet[row_i, ][["integration.case"]]
        path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
        
        path.to.save.html <- file.path(output.dir, sprintf("from_%s", output.index), integration.case, sprintf("%s_%s", sample1, sample2))
        output.file.name <- sprintf("%s_vs_%s.part1.html", sample1, sample2)
        
        path.to.save.DGE.output <- file.path(path.to.13.output, sprintf("from_%s", output.index), integration.case, sprintf("%s_%s", sample1, sample2), "part1")
        
        input.params <- list(
          sample1 = sample1,
          sample2 = sample2,
          path.to.s.obj = path.to.s.obj,
          path.to.save.output = path.to.save.DGE.output
        )
      } else if (output.index == "10_output"){
        integration.case <- input.samplesheet[row_i, ][["integration.case"]]
        regression.mode <- input.samplesheet[row_i, ][["regression.mode"]]
        filter.mode <- input.samplesheet[row_i, ][["filter.mode"]]
        path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
        
        path.to.save.html <- file.path(output.dir, sprintf("from_%s", output.index), integration.case, regression.mode, filter.mode, sprintf("%s_%s", sample1, sample2))
        output.file.name <- sprintf("%s_vs_%s.part1.html", sample1, sample2)
        
        path.to.save.DGE.output <- file.path(path.to.13.output, sprintf("from_%s", output.index), integration.case, regression.mode, filter.mode, sprintf("%s_%s", sample1, sample2), "part1")
        
        input.params <- list(
          sample1 = sample1,
          sample2 = sample2,
          path.to.s.obj = path.to.s.obj,
          path.to.save.output = path.to.save.DGE.output
        )
      }
      rmarkdown::render(path.to.rmd,
                        params = input.params,
                        output_file = output.file.name,
                        output_dir = path.to.save.html)    
    }
  }
}

