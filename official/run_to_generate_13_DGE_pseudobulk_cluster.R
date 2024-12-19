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
path.to.rmd <- file.path(path.to.main.src, src.dir, "13_DGE_analysis_pseudobulk_clusters.Rmd")
output.dir <- file.path(path.to.save.html, "13_output")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

samplesheets <- list(
  `03_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_03_output_simplified.xlsx")),
  `10_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_10_output_simplified.xlsx")),
  `12_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_12_output_simplified.xlsx")) %>%
    subset(output_index == "12_output"),
  `12_output_remove_BCR_TCR` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_12_output_simplified.xlsx")) %>%
    subset(output_index == "12_output_remove_BCR_TCR"),
  `17_output` = readxl::read_excel(file.path(path.to.save.samplesheet, "SampleSheet_17_output_simplified.xlsx"))
)

#####----------------------------------------------------------------------------#####

comparison.samplesheet <- read.csv(file.path(path.to.main.src, src.dir, "sample_comparision_list.csv"))

for (output.index in names(samplesheets)){
  input.samplesheet <- samplesheets[[output.index]]
  for (row_i in seq(1, nrow(input.samplesheet))){
    integration.case <- input.samplesheet[row_i, ][["integration.case"]]
    path.to.s.obj <- input.samplesheet[row_i, ][["path"]]
    s.obj <- readRDS(path.to.s.obj)
    meta.data <- s.obj@meta.data %>%
      rowwise() %>%
      mutate(name_ht = sprintf("%s_%s", name, HTO_classification))
    countdf <- table(meta.data$name_ht, meta.data$cca.cluster.0.5) %>% as.data.frame() %>%
      subset(grepl("Hashtag", Var1) == TRUE) %>%
      pivot_wider(names_from = "Var1", values_from = "Freq") %>% column_to_rownames("Var2")
    
    for (j in seq(1, nrow(comparison.samplesheet))){
      sample1 <- comparison.samplesheet[j, ][["sample1"]]
      sample2 <- comparison.samplesheet[j, ][["sample2"]]
      tmp.countdf <- countdf[, to_vec(for (item in colnames(countdf)) if(str_split(item, "_Hashtag")[[1]][[1]] == sample1 | str_split(item, "_Hashtag")[[1]][[1]] == sample2) item)]
      condition <- function(row) {
        all(row > 1)
      }
      tmp.countdf$check <- apply(tmp.countdf, 1, condition)
      tmp.countdf$count.sample1 <- unlist(lapply(
        row.names(tmp.countdf), function(x){
          to_vec(for (item in names(tmp.countdf[x, ])) if(sample1 == str_split(item, "_Hashtag")[[1]][[1]] & tmp.countdf[x, ][[item]] != 0) item) %>% length()
        }
      )) 
      tmp.countdf$count.sample2 <- unlist(lapply(
        row.names(tmp.countdf), function(x){
          to_vec(for (item in names(tmp.countdf[x, ])) if(sample2 == str_split(item, "_Hashtag")[[1]][[1]] & tmp.countdf[x, ][[item]] != 0) item) %>% length()
        }
      )) 
      available.clusters <- row.names(subset(tmp.countdf, tmp.countdf$count.sample1 > 1 & tmp.countdf$count.sample2 > 1 & tmp.countdf$check == TRUE))
      for (cluster.id in available.clusters){
        if (output.index == "03_output"){
          path.to.save.html <- file.path(output.dir, 
                                         sprintf("from_%s", output.index), 
                                         integration.case, 
                                         sprintf("%s_%s", sample1, sample2),
                                         sprintf("cluster_%s", cluster.id))
          output.file.name <- sprintf("%s_vs_%s.part2.html", sample1, sample2)
          
          path.to.save.DGE.output <- file.path(path.to.13.output, 
                                               sprintf("from_%s", output.index), 
                                               integration.case, 
                                               sprintf("%s_%s", sample1, sample2), 
                                               "part2",
                                               sprintf("cluster_%s", cluster.id))
          input.params <- list(
            sample1 = sample1,
            sample2 = sample2,
            path.to.s.obj = path.to.s.obj,
            path.to.save.output = path.to.save.DGE.output,
            cluster.id = cluster.id
          )
        } else if (output.index == "10_output"){
          regression.mode <- input.samplesheet[row_i, ][["regression.mode"]]
          filter.mode <- input.samplesheet[row_i, ][["filter.mode"]]
          
          path.to.save.html <- file.path(output.dir, sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sprintf("%s_%s", sample1, sample2),
                                         sprintf("cluster_%s", cluster.id))
          output.file.name <- sprintf("%s_vs_%s.part2.html", sample1, sample2)
          
          path.to.save.DGE.output <- file.path(path.to.13.output, 
                                               sprintf("from_%s", output.index), 
                                               integration.case, 
                                               regression.mode, 
                                               filter.mode, 
                                               sprintf("%s_%s", sample1, sample2), 
                                               "part2",
                                               sprintf("cluster_%s", cluster.id))
          
          input.params <- list(
            sample1 = sample1,
            sample2 = sample2,
            path.to.s.obj = path.to.s.obj,
            path.to.save.output = path.to.save.DGE.output,
            cluster.id = cluster.id
          )
        } else if (output.index %in% c("12_output", "12_output_remove_BCR_TCR", "17_output")){
          regression.mode <- input.samplesheet[row_i, ][["regression.mode"]]
          filter.mode <- input.samplesheet[row_i, ][["filter.mode"]]
          sub.cluster.id <- input.samplesheet[row_i, ][["sub.cluster.id"]]
          
          path.to.save.html <- file.path(output.dir, sprintf("from_%s", output.index), 
                                         integration.case, 
                                         regression.mode, 
                                         filter.mode, 
                                         sub.cluster.id,
                                         sprintf("%s_%s", sample1, sample2),
                                         sprintf("cluster_%s", cluster.id))
          output.file.name <- sprintf("%s_vs_%s.part2.html", sample1, sample2)
          
          path.to.save.DGE.output <- file.path(path.to.13.output, 
                                               sprintf("from_%s", output.index), 
                                               integration.case, 
                                               regression.mode, 
                                               filter.mode, 
                                               sub.cluster.id,
                                               sprintf("%s_%s", sample1, sample2), 
                                               "part2",
                                               sprintf("cluster_%s", cluster.id))
          input.params <- list(
            sample1 = sample1,
            sample2 = sample2,
            path.to.s.obj = path.to.s.obj,
            path.to.save.output = path.to.save.DGE.output,
            cluster.id = cluster.id
          )
        }
        
        dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
        for (param.name in names(input.params)){
          print(sprintf("Param %s: %s", param.name, input.params[[param.name]]))
        }
        if (file.exists(file.path(path.to.save.html, output.file.name)) ==  FALSE){
          print(sprintf("Generating html file output %s", file.path(path.to.save.html, output.file.name)))
          rmarkdown::render(path.to.rmd,
                            params = input.params,
                            output_file = output.file.name,
                            output_dir = path.to.save.html)          
        } else {
          print(sprintf("File %s exists", file.path(path.to.save.html, output.file.name)))
        }
      }      
    }
  }
}

