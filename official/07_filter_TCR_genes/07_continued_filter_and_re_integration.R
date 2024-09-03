gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))

options(future.globals.maxSize = 10000 * 1024^2)

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
# filter.mode <- "TCR_genes" # filter.mode could be "TCR_genes", "nCounts" or "nCount_and_TCRgenes"
# integration.case <- "remove_d4_LPS"
# regression.mode <- "CC_differences"
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

#####----------------------------------------------------------------------#####
##### SOME CONFIGURATIONS
#####----------------------------------------------------------------------#####
integration.config <- list(
  all.samples = c("adult_GF",
                  "d4_GF",
                  "adult_SPF",
                  "d4_LPS",
                  "d10_SPF",
                  "d4_SPF",
                  "d15_SPF",
                  "d7_GF",
                  "d20_SPF",
                  "d7_SPF",
                  "SC12",
                  "SC11",
                  "SC5"),
  remove_d4_LPS_SC5 = c("adult_GF",
                        "d4_GF",
                        "adult_SPF",
                        "d10_SPF",
                        "d4_SPF",
                        "d15_SPF",
                        "d7_GF",
                        "d20_SPF",
                        "d7_SPF",
                        "SC12",
                        "SC11"),
  remove_d4_LPS = c("adult_GF",
                    "d4_GF",
                    "adult_SPF",
                    "d10_SPF",
                    "d4_SPF",
                    "d15_SPF",
                    "d7_GF",
                    "d20_SPF",
                    "d7_SPF",
                    "SC12",
                    "SC11",
                    "SC5")
)

cell.cycle.features <- list(
  CC_differences = c("CC.Difference"),
  S_G2M_G1_scores = c("S.Score", "G2M.Score", "G1.Score")
)

TR_genes_patterns <- c("Trav", "Traj", "Trac", "Trbv", "Trbd", "Trbj", "Trbc",
                       "Trgv", "Trgj", "Trgc", "Trdv", "Trdc", "Trdj") 
BR_genes_patterns <- c("Ighv", "Ighd", "Ighj", "Ighc", "Igkv",
                         "Igkj", "Igkc", "Iglv", "Iglj", "Iglc")
#####----------------------------------------------------------------------#####
##### Define functions
#####----------------------------------------------------------------------#####
run_07 <- function(outdir, 
                   PROJECT, 
                   integration.case, 
                   regression.mode, 
                   filter.mode, 
                   integration.config, 
                   cell.cycle.features, 
                   TR_genes_patterns){
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  path.to.07.output <- file.path(path.to.main.output, "07_output", integration.case, regression.mode, filter.mode)
  dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)
  
  integrate.samples <- integration.config[[integration.case]]
  
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  regressOut_mode <- NULL
  features_to_regressOut <- NULL
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt", cell.cycle.features[[regression.mode]])
  cluster.resolution <- 0.5
  
  if (file.exists(file.path(path.to.07.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.csv")) == FALSE){
    if (file.exists(file.path(path.to.07.output, sprintf("raw_merge_dataset_%s.%s.%s.csv", integration.case, regression.mode, filter.mode))) == FALSE){
      data.list <- list()
      for (i in seq(length(integrate.samples))){
        sample.id <- integrate.samples[[i]]
        path.to.05.output <- file.path(path.to.main.output, "05_output", sample.id, regression.mode)
        print(sprintf("reading in sample %s", sample.id))
        tmp.s.obj <- readRDS(file.path(path.to.05.output, sprintf("%s.cellcycle_reg.rds", sample.id)))
        
        #####--------------------------------------------------------------------#####
        ##### DO FILTERING STEPS HERE
        #####--------------------------------------------------------------------#####
        
        ##### remove TCR genes
        all.genes <- row.names(tmp.s.obj)
        TCRgenes.to.exclude <- unlist(lapply(all.genes, function(x){
          if (substr(x, 1, 4) %in% TR_genes_patterns){
            return(x)
          } else {
            return(NA)
          }
        }))
        TCRgenes.to.exclude <- subset(TCRgenes.to.exclude, is.na(TCRgenes.to.exclude) == FALSE)
        
        ##### remove BCR genes 
        BCRgenes.to.exclude <- unlist(lapply(all.genes, function(x){
          if (substr(x, 1, 4) %in% BR_genes_patterns){
            return(x)
          } else {
            return(NA)
          }
        }))
        BCRgenes.to.exclude <- subset(BCRgenes.to.exclude, is.na(BCRgenes.to.exclude) == FALSE)

        ##### remove doublet cells by counts
        nCount <- tmp.s.obj@meta.data$nCount_RNA
        tmp.metadata <- tmp.s.obj@meta.data %>% rownames_to_column("barcode") %>%
          mutate(complexity = log10(nFeature_RNA) / log10(nCount_RNA) )
        all.qt <- quantile(nCount) %>% as.list()
        iqr <- all.qt$`75%` - all.qt$`25%`
        upper.cut <- all.qt$`75%` + 2.5 * iqr
        selected.cells <- subset(tmp.metadata, tmp.metadata$nCount_RNA <= upper.cut & tmp.metadata$complexity >= 0.8)
        
        if (filter.mode == "TCR_genes"){
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, TCRgenes.to.exclude))
        } else if (filter.mode == "nCount"){
          tmp.s.obj <- subset(tmp.s.obj, cells = selected.cells$barcode)
        } else if (filter.mode == "nCount_and_TCRgenes"){
          tmp.s.obj <- subset(tmp.s.obj, cells = selected.cells$barcode)
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, TCRgenes.to.exclude))
        } else if (filter.mode == "BCR_genes"){
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, BCRgenes.to.exclude))
        } else if (filter.mode == "BCR_TCR_genes"){
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, TCRgenes.to.exclude))
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, BCRgenes.to.exclude))
        } else if (filter.mode == "nCount_and_BCR_TCRgenes"){
          tmp.s.obj <- subset(tmp.s.obj, cells = selected.cells$barcode)
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, TCRgenes.to.exclude))
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, BCRgenes.to.exclude))
        } else if (filter.mode == "nCount_and_BCRgenes"){
          tmp.s.obj <- subset(tmp.s.obj, cells = selected.cells$barcode)
          tmp.s.obj <- subset(tmp.s.obj, features = setdiff(all.genes, BCRgenes.to.exclude))
        }
        data.list[[i]] <- tmp.s.obj
      }
      
      s.obj <- merge(data.list[[1]], data.list[2: length(integrate.samples)])
      saveRDS(s.obj, file.path(path.to.07.output, sprintf("raw_merge_dataset_%s.%s.%s.rds", integration.case, regression.mode, filter.mode)))
      write.csv(data.frame(status = c("finished saving file")), 
                file.path(path.to.07.output, sprintf("raw_merge_dataset_%s.%s.%s.csv", integration.case, regression.mode, filter.mode)))
    } else {
      print("reading data in to s.obj ...")
      s.obj <- readRDS(file.path(path.to.07.output, sprintf("raw_merge_dataset_%s.%s.%s.rds", integration.case, regression.mode, filter.mode)))
      print("finished reading in saved data ...")
    }
    
    DefaultAssay(s.obj) <- "RNA"
    print("Running JoinLayers ...")
    s.obj <- JoinLayers(s.obj)
    print("Finished JoinLayers,...")
    print("Start integration")
    s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                         save.RDS.s8 = TRUE,
                                                         path.to.output = path.to.07.output,
                                                         use.sctransform = TRUE,
                                                         num.PCA = num.PCA,
                                                         num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                         num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                         cluster.resolution = cluster.resolution,
                                                         vars.to.regress = vars.to.regress)        
    write.csv(data.frame(status = c("finished saving...")), 
              file.path(path.to.07.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.csv"))
  } else {
    print(sprintf("integrated done. file at %s", file.path(path.to.07.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.rds")))
  }
}

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####
for (filter.mode in c("TCR_genes",
                      "nCount",
                      "nCount_and_TCRgenes",
                      "BCR_genes",
                      "nCount_and_BCRgenes",
                      "nCount_and_BCR_TCRgenes",
                      "BCR_TCR_genes")){
  for (regression.mode in names(cell.cycle.features)){
    for (integration.case in names(integration.config)){
      print(sprintf("Working on %s - %s - %s", filter.mode, regression.mode, integration.case))
      run_07(outdir = outdir, 
             PROJECT = PROJECT, 
             integration.case = integration.case, 
             regression.mode = regression.mode, 
             filter.mode = filter.mode, 
             integration.config = integration.config, 
             cell.cycle.features = cell.cycle.features, 
             TR_genes_patterns = TR_genes_patterns)
    }
  }
}
