#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs", "13_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "13_DGE"
path.to.rmd <- file.path(path.to.main.src, src.dir, "13_DGE_analysis.Rmd")
path.to.pseudobulk.rmd <- file.path(path.to.main.src, src.dir, "13_DGE_analysis_pseudobulk_clusters.Rmd")

path.to.13.output <- file.path(path.to.save.html, "13_output")
dir.create(path.to.13.output, showWarnings = FALSE, recursive = TRUE)

samplesheet <- readxl::read_excel(file.path(path.to.main.src, "SampleSheet_for_DGE_CellChat_Monocle_RNAvelocity.xlsx"))
samplelist <- read.csv(file.path(path.to.main.src, src.dir, "sample_comparision_list.csv"))

for (i in seq(1, nrow(samplesheet))){
  for (j in seq(1, nrow(samplelist))){
    input.info <- c()
    input.info.cols <- c("integration.case", 
                         "regression.mode", 
                         "filter.mode", 
                         "sub.cluster.id")
    
    for (c in input.info.cols){
      if (is.na(samplesheet[i, ][[c]]) == FALSE){
        input.info <- c(input.info, samplesheet[i, ][[c]])
      }
    }
    
    sample1 <- samplelist[j, ]$sample1
    sample2 <- samplelist[j, ]$sample2
    path.to.s.obj <- samplesheet[i, ]$path
    
    path.to.save.output <- file.path(
      path.to.13.output, 
      sprintf("from_%s", samplesheet[i, ]$output_index),
      paste0(input.info, collapse = "/"),
      sprintf("%s_vs_%s", sample1, sample2)
    )
    
    output_file <- sprintf("%s_vs_%s.part1.html", sample1, sample2)
    output_dir <- file.path(
      path.to.save.html, 
      sprintf("from_%s", samplesheet[i, ]$output_index),
      paste0(input.info, collapse = "/"),
      sprintf("%s_vs_%s", sample1, sample2)
    )
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ##### Generate single cell DGE
    if (file.exists(file.path(output_dir, output_file)) == FALSE){
      print(sprintf("Generate new html file %s...", file.path(output_dir, output_file) ))
      rmarkdown::render(input = path.to.rmd,
                        output_file = output_file,
                        output_dir = output_dir,
                        params = list(
                          sample1 = sample1,
                          sample2 = sample2,
                          path.to.save.output = path.to.save.output,
                          path.to.s.obj = path.to.s.obj))
    } else {
      print(sprintf("html file %s exists...", file.path(output_dir, output_file) ))
    }
    
    ##### Generate pseudobulk DGE
    s.obj <- readRDS(path.to.s.obj)
    num.clusters <- length(unique(s.obj$cca.cluster.0.5))
    print(sprintf("Number of clusters to perform pseudobulk DGE: %s", num.clusters))
    
    meta.data <- s.obj@meta.data %>%
      rowwise() %>%
      mutate(name_ht = sprintf("%s_%s", name, HTO_classification))
    countdf <- table(meta.data$name_ht, meta.data$cca.cluster.0.5) %>% as.data.frame() %>%
      subset(grepl("Hashtag", Var1) == TRUE) %>%
      pivot_wider(names_from = "Var1", values_from = "Freq") %>% column_to_rownames("Var2")
    
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
    
    for (k in available.clusters){
      output_dir2 <- file.path(
        path.to.save.html, 
        sprintf("from_%s", samplesheet[i, ]$output_index),
        paste0(input.info, collapse = "/"),
        sprintf("%s_vs_%s", sample1, sample2),
        sprintf("cluster_%s", k)
      )
      path.to.save.pseudobulk.output <- file.path(
        path.to.13.output, 
        sprintf("from_%s", samplesheet[i, ]$output_index),
        paste0(input.info, collapse = "/"),
        sprintf("%s_vs_%s", sample1, sample2),
        "part2",
        sprintf("cluster_%s", k)
      )
      dir.create(output_dir2, showWarnings = FALSE, recursive = TRUE)
      output_file2 <- sprintf("%s_vs_%s.part2.html", sample1, sample2)
      
      if (file.exists(file.path(output_dir2, output_file2)) == FALSE){
        print(sprintf("Generate new html file %s...", file.path(output_dir, output_file) ))
        rmarkdown::render(input = path.to.pseudobulk.rmd,
                          output_file = output_file2,
                          output_dir = output_dir2,
                          params = list(
                            sample1 = sample1,
                            sample2 = sample2,
                            path.to.save.output = path.to.save.pseudobulk.output,
                            path.to.s.obj = path.to.s.obj,
                            cluster.id = k))
      } else {
        print(sprintf("html file %s exists...", file.path(output_dir, output_file) ))
      }
    }
  }
}

