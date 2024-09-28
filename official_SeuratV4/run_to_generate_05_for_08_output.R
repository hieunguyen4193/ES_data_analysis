#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 04 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
s.obj.name <- "all_sobj.integrated.rds"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.08.output <- file.path(path.to.main.output, "data_analysis", "08_output")
path.to.05.output <- file.path(path.to.main.output, "data_analysis", "05_output")
path.to.04.output <- file.path(path.to.main.output, "data_analysis", "04_output")

path.to.save.html <- file.path(path.to.main.output, "html_outputs", "05_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

all.path.to.s.obj <- list(
  `04_output` = list(
    Myeloid_Basophils = file.path(path.to.04.output, "all_sobj.integrated.rds", "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", "Myeloid_Basophils")),
    T_cells = file.path(path.to.04.output, "all_sobj.integrated.rds", "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", "T_cells")),
    B_cells = file.path(path.to.04.output, "all_sobj.integrated.rds", "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", "B_cells")),
    full = file.path(path.to.04.output, "d4_samples_only", "all_sobj.integrated_d4_LPS_d4_SPF_only.rds")
  ),
  `08_output` = list(
    Myeloid_Basophils = list(
      ambientRNAonly = file.path(path.to.08.output, "Myeloid_Basophils", "Myeloid_Basophils.remove_contam_T_B_cells.ambientRNAonly.rds"),
      all_cells = file.path(path.to.08.output, "Myeloid_Basophils", "Myeloid_Basophils.remove_contam_T_B_cells.rds")
    )
  )
)

#####----------------------------------------------------------------------#####
##### Run to generate DGE in 05_DGE_analysis.Rmd
#####----------------------------------------------------------------------#####
for (output.index in names(all.path.to.s.obj)){
  for (celltype in names(all.path.to.s.obj[[output.index]])){
    for (sub.cluster.idx in names(all.path.to.s.obj[[output.index]][[celltype]])){
      path.to.s.obj <- all.path.to.s.obj[[output.index]][[celltype]][[sub.cluster.idx]]
      path.to.rmd <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4/05_DGE_analysis.Rmd"
      output_file <- sprintf("05_DGE_%s_vs_%s.html", sample1, sample2)
      path.to.save.output <- file.path(path.to.05.output, output.index, celltype, sub.cluster.idx)
      output_dir <- file.path(path.to.save.html, output.index, celltype, sub.cluster.idx)
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      if (file.exists(file.path(output_dir, output_file)) == FALSE){
        rmarkdown::render(input = path.to.rmd, 
                          output_file = output_file,
                          output_dir = output_dir,
                          params = list(sample1= sample1,
                                        sample2= sample2,
                                        path.to.s.obj= path.to.s.obj,
                                        path.to.save.output= path.to.save.output))    
      }
    }
  }
}

#####----------------------------------------------------------------------#####
##### Run to generate DGE in 05_DGE_analysis.pseudoBulk_cluster.Rmd
#####----------------------------------------------------------------------#####
for (output.index in names(all.path.to.s.obj)){
  
  for (celltype in names(all.path.to.s.obj[[output.index]])){
    for (sub.cluster.idx in names(all.path.to.s.obj[[output.index]][[celltype]])){
      path.to.s.obj <- all.path.to.s.obj[[output.index]][[celltype]][[sub.cluster.idx]]
      s.obj <- readRDS(path.to.s.obj)
      meta.data <- s.obj@meta.data %>%
        rowwise() %>%
        mutate(name_ht = sprintf("%s_%s", name, HTO_classification))
      countdf <- table(meta.data$name_ht, meta.data$seurat_clusters) %>% as.data.frame() %>%
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
      for (cluster.id in available.clusters){
        path.to.rmd <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4/05_DGE_analysis.pseudoBulk_cluster.Rmd"
        output_file <- sprintf("05_DGE_%s_vs_%s.cluster_%s.html", sample1, sample2, cluster.id)
        path.to.save.output <- file.path(path.to.05.output, output.index, celltype, sub.cluster.idx, cluster.id)
        output_dir <- file.path(path.to.save.html, output.index, celltype, sub.cluster.idx, cluster.id)
        dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        if (file.exists(file.path(output_dir, output_file)) == FALSE){
          rmarkdown::render(input = path.to.rmd, 
                            output_file = output_file,
                            output_dir = output_dir,
                            params = list(sample1= sample1,
                                          sample2= sample2,
                                          path.to.s.obj= path.to.s.obj,
                                          path.to.save.output= path.to.save.output,
                                          cluster.id = cluster.id))  
        }
      }
    }
  }
  
}