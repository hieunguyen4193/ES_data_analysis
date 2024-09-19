#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

##### install CellChat
if ("CellChat" %in% row.names(installed.packages()) == FALSE){
  devtools::install_github("immunogenomics/presto", upgrade = "never")
  devtools::install_github("jinworks/CellChat", upgrade = "never")
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", repos = NULL, type = "source")
} 

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.14.output <- file.path(outdir, PROJECT, "data_analysis", "14_output")

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path.to.save.samplesheet <- file.path(path.to.main.src, "sampleSheets_for_DGE_and_CellChat")

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "14_CellChat"
path.to.rmd <- file.path(path.to.main.src, src.dir, "14_CellChat_general_analysis.Rmd")
output.dir <- file.path(path.to.save.html, "14_output")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
  `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output_simplified.xlsx")),
  `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output_simplified.xlsx")),
  `12_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_12_output_simplified.xlsx"))
)

all.samples <- c("adult_GF",
                 "adult_SPF",
                 "d10_SPF",
                 "d15_SPF",
                 "d20_SPF",
                 "d4_GF",
                 "d4_SPF",
                 "d7_GF",
                 "d7_SPF",
                 "SC11",
                 "SC12")

filter10cells <- "Filter10"
for (sample.id in all.samples){
  for (output.index in names(samplesheets)){
    input.samplesheet <- samplesheets[[output.index]]
    for (row_i in seq(1, nrow(input.samplesheet))){
      if (output.index == "03_output"){
        integration.case <- input.samplesheet[row_i, ][["integration.case"]]
        path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
        
        path.to.save.html <- file.path(output.dir, 
                                       sprintf("from_%s", output.index), 
                                       integration.case, 
                                       sample.id)
        output.file.name <- sprintf("%s.CellChat_single_sample.html", sample.id)
        
        path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                  sprintf("from_%s", output.index), 
                                                  integration.case, 
                                                  sample.id, 
                                                  "part1")
        
        input.params <- list(
          sample.id = sample.id,
          filter10cells = filter10cells,
          path.to.s.obj = path.to.s.obj,
          path.to.save.output = path.to.save.CellChat.output
        )
      } else if (output.index == "10_output"){
        integration.case <- input.samplesheet[row_i, ][["integration.case"]]
        regression.mode <- input.samplesheet[row_i, ][["regression.mode"]]
        filter.mode <- input.samplesheet[row_i, ][["filter.mode"]]
        path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
        
        path.to.save.html <- file.path(output.dir, 
                                       sprintf("from_%s", output.index), 
                                       integration.case, 
                                       regression.mode, 
                                       filter.mode, 
                                       sample.id)
        output.file.name <- sprintf("%s.CellChat_single_sample.html", sample.id)
        
        path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                  sprintf("from_%s", output.index), 
                                                  integration.case, 
                                                  regression.mode, 
                                                  filter.mode, 
                                                  sample.id, 
                                                  "part1")
        
        input.params <- list(
          sample.id = sample.id,
          filter10cells = filter10cells,
          path.to.s.obj = path.to.s.obj,
          path.to.save.output = path.to.save.CellChat.output
        )
      } else if (output.index == "12_output"){
        integration.case <- input.samplesheet[row_i, ][["integration.case"]]
        regression.mode <- input.samplesheet[row_i, ][["regression.mode"]]
        filter.mode <- input.samplesheet[row_i, ][["filter.mode"]]
        sub.cluster.id <- input.samplesheet[row_i, ][["sub.cluster.id"]]
        
        path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
        
        path.to.save.html <- file.path(output.dir, 
                                       sprintf("from_%s", output.index), 
                                       integration.case, 
                                       regression.mode, 
                                       filter.mode, 
                                       sub.cluster.id,
                                       sample.id)
        output.file.name <- sprintf("%s.%s.CellChat_single_sample.html", sample.id, sub.cluster.id)
        
        path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                  sprintf("from_%s", output.index), 
                                                  integration.case, 
                                                  regression.mode, 
                                                  filter.mode, 
                                                  sub.cluster.id,
                                                  sample.id, 
                                                  "part1")
        
        input.params <- list(
          sample.id = sample.id,
          filter10cells = filter10cells,
          path.to.s.obj = path.to.s.obj,
          path.to.save.output = path.to.save.CellChat.output
        )
      }
      for (param.name in names(input.params)){
        print(sprintf("Param %s: %s", param.name, input.params[[param.name]]))
      }
      if (file.exists(file.path(path.to.save.html, output.file.name)) == FALSE){
        rmarkdown::render(path.to.rmd,
                          params = input.params,
                          output_file = output.file.name,
                          output_dir = path.to.save.html)            
      } else {
        print(sprintf("File html %s already exists, not generating new file.", file.path(path.to.save.html, output.file.name)))
      }
    }
  }
}
