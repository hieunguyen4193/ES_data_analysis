gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))

options(future.globals.maxSize = 10000 * 1024^2)

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"
chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25
num.PC.used.in.UMAP <- 25
my_random_seed <- 42

#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####
# PROJECT <- params$PROJECT
# sample.id <- params$sample.id
# outdir <- params$outdir

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV5"

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
                    "SC5"),
  remove_d4_LPS_SC5_SC11 = c("adult_GF",
                        "d4_GF",
                        "adult_SPF",
                        "d10_SPF",
                        "d4_SPF",
                        "d15_SPF",
                        "d7_GF",
                        "d20_SPF",
                        "d7_SPF",
                        "SC12")
)

for (integration.case in names(integration.config)){
  
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.02.output <- file.path(path.to.main.output, "02_output")
  path.to.03.output <- file.path(path.to.main.output, "03_output", integration.case)
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
    print(sprintf("working on sample %s", integration.case))
    integrate.samples <- integration.config[[integration.case]]
    if (file.exists(file.path(path.to.03.output, sprintf("raw_merge_dataset_%s.rds", integration.case))) == FALSE){
      data.list <- list()
      for (i in seq(length(integrate.samples))){
        sample.id <- integrate.samples[[i]]
        print(sprintf("reading in sample %s", sample.id))
        data.list[[i]] <- readRDS(file.path(path.to.02.output, sample.id, sprintf("GEX_sample_%s_seurat_object.rds", sample.id)))
      }
      s.obj <- merge(data.list[[1]], data.list[2: length(integrate.samples)])
      saveRDS(s.obj, file.path(path.to.03.output, sprintf("raw_merge_dataset_%s.rds", integration.case)))  
    } else {
      s.obj <- readRDS(file.path(path.to.03.output, sprintf("raw_merge_dataset_%s.rds", integration.case)))
    }
    
    num.PCA <- 25
    num.PC.used.in.UMAP <- 25
    num.PC.used.in.Clustering <- 25
    regressOut_mode <- NULL
    features_to_regressOut <- NULL
    use.sctransform <- TRUE
    vars.to.regress <- c("percent.mt")
    cluster.resolution <- 0.5
    
    DefaultAssay(s.obj) <- "RNA"
    s.obj <- JoinLayers(s.obj)
    s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                         save.RDS.s8 = TRUE,
                                                         path.to.output = path.to.03.output,
                                                         use.sctransform = TRUE,
                                                         num.PCA = num.PCA,
                                                         num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                         num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                         cluster.resolution = cluster.resolution,
                                                         vars.to.regress = vars.to.regress)
  }
}
