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

# get this path to get Cellchat data object for a single sample
path.to.14.output.single <- file.path(outdir, PROJECT, "data_analysis", "14_output")

path.to.14.output <- file.path(outdir, PROJECT, "data_analysis", "14_output", "compare_2_samples")

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path.to.save.samplesheet <- file.path(path.to.main.src, "sampleSheets_for_DGE_and_CellChat")

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "14_CellChat"
path.to.rmd <- file.path(path.to.main.src, src.dir, "14_CellChat_diff_analysis_2_samples.Rmd")
output.dir <- file.path(path.to.save.html, "14_output")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
  `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output_simplified.xlsx")),
  `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output_simplified.xlsx")),
  `12_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_12_output_simplified.xlsx")) %>%
    subset(output_index == "12_output"),
  `12_output_remove_BCR_TCR` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_12_output_simplified.xlsx")) %>%
    subset(output_index == "12_output_remove_BCR_TCR")
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
comparison.samplesheet <- read.csv(file.path(path.to.main.src, "13_DGE", "sample_comparision_list.csv"))

rerun <- FALSE
for (i in seq(1, nrow(comparison.samplesheet))){
  sample1 <- comparison.samplesheet[i, "sample1"]
  sample2 <- comparison.samplesheet[i, "sample2"]
  if (sample1 %in% c("d10_SPF", "adult_SPF") & sample2 == "SC12"){
  print(sprintf("skip this case, when sample1 = %s and sample2 = %s", sample1, sample2))
  } else {
    for (output.index in names(samplesheets)){
      input.samplesheet <- samplesheets[[output.index]]
      for (row_i in seq(1, nrow(input.samplesheet))){
        if (output.index == "03_output"){
          integration.case <- input.samplesheet[row_i, ][["integration.case"]]
          path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
          
          path.to.save.html <- file.path(output.dir, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         sprintf("%s_vs_%s", sample1, sample2))
          output.file.name <- sprintf("CellChat_%s_vs_%s.html", sample1, sample2)
          
          path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                    sprintf("from_%s", output.index), 
                                                    integration.case, 
                                                    sprintf("%s_vs_%s", sample1, sample2))
          path.to.cellchat1 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         sample1,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample1))
          
          path.to.cellchat2 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         sample2,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample2))
          
          input.params <- list(
            path.to.s.obj = path.to.s.obj,
            sample1 = sample1,
            path.to.cellchat1 = path.to.cellchat1,
            sample2 = sample2,
            path.to.cellchat2 = path.to.cellchat2,
            path.to.save.output = path.to.save.CellChat.output,
            filter10cells = filter10cells
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
                                         sprintf("%s_vs_%s", sample1, sample2))
          output.file.name <- sprintf("CellChat_%s_vs_%s.html", sample1, sample2)
          
          path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                    sprintf("from_%s", output.index), 
                                                    integration.case, 
                                                    regression.mode, 
                                                    filter.mode, 
                                                    sprintf("%s_vs_%s", sample1, sample2))
          path.to.cellchat1 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sample1,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample1))
          
          path.to.cellchat2 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sample2,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample2))
          input.params <- list(
            path.to.s.obj = path.to.s.obj,
            sample1 = sample1,
            path.to.cellchat1 = path.to.cellchat1,
            sample2 = sample2,
            path.to.cellchat2 = path.to.cellchat2,
            path.to.save.output = path.to.save.CellChat.output,
            filter10cells = filter10cells
          )
        } else if (output.index %in% c("12_output", "12_output_remove_BCR_TCR")){
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
                                         sprintf("%s_vs_%s", sample1, sample2))
          output.file.name <- sprintf("CellChat_%s_%s_vs_%s.html", sub.cluster.id, sample1, sample2)
          
          path.to.save.CellChat.output <- file.path(path.to.14.output, 
                                                    sprintf("from_%s", output.index), 
                                                    integration.case, 
                                                    regression.mode, 
                                                    filter.mode, 
                                                    sub.cluster.id,
                                                    sprintf("%s_vs_%s", sample1, sample2))
          path.to.cellchat1 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sub.cluster.id,
                                         sample1,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample1))
          
          path.to.cellchat2 <- file.path(path.to.14.output.single, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sub.cluster.id,
                                         sample2,
                                         "part1", 
                                         sprintf("CellChat_object.%s.Filter10.rds", sample2))
          input.params <- list(
            path.to.s.obj = path.to.s.obj,
            sample1 = sample1,
            path.to.cellchat1 = path.to.cellchat1,
            sample2 = sample2,
            path.to.cellchat2 = path.to.cellchat2,
            path.to.save.output = path.to.save.CellChat.output,
            filter10cells = filter10cells
          )
        }
        for (param.name in names(input.params)){
          print(sprintf("Param %s: %s", param.name, input.params[[param.name]]))
        }
        if (file.exists(file.path(path.to.save.html, output.file.name)) == FALSE | rerun == TRUE){
          print(sprintf("Generate the file %s ...", output.file.name))
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
}


