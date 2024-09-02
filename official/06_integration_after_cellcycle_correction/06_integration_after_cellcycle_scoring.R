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
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

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

cell.cycle.features <- list(
  CC_differences = c("CC.Difference"),
  S_G2M_G1_scores = c("S.Score", "G2M.Score", "G1.Score")
)

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.merge.data <- file.path(path.to.main.output, "06_output", "merge_data")
dir.create(path.to.merge.data, showWarnings = FALSE, recursive = TRUE)

for (regression.mode in c("CC_differences", "S_G2M_G1_scores")){
  for (integration.case in names(integration.config)){
    
    path.to.06.output <- file.path(path.to.main.output, "06_output", integration.case, regression.mode)
    dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
    
    if (file.exists(file.path(path.to.06.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
      print(sprintf("working on sample %s", integration.case))
      integrate.samples <- integration.config[[integration.case]]
      
      ##### merging the data before integration
      if (file.exists(file.path(path.to.merge.data, sprintf("raw_merge_dataset_%s.rds", integration.case))) == FALSE){
        print("Merge data does not exists, generate new merged data...")
        data.list <- list()
        for (i in seq(length(integrate.samples))){
          sample.id <- integrate.samples[[i]]
          path.to.05.output <- file.path(path.to.main.output, "05_output", sample.id, regression.mode)
          print(sprintf("reading in sample %s", sample.id))
          data.list[[i]] <- readRDS(file.path(path.to.05.output, sprintf("%s.cellcycle_reg.rds", sample.id)))
        }
        s.obj <- merge(data.list[[1]], data.list[2: length(integrate.samples)])
        saveRDS(s.obj, file.path(path.to.merge.data, sprintf("raw_merge_dataset_%s.rds", integration.case)))  
        write.csv(data.frame(status = c("finished saving data")), file.path(path.to.06.output, "finished_saving_data.csv"))
      } else {
        print("reading in raw merge data ...")
        s.obj <- readRDS(file.path(path.to.merge.data, sprintf("raw_merge_dataset_%s.rds", integration.case)))
      }
      
      num.PCA <- 25
      num.PC.used.in.UMAP <- 25
      num.PC.used.in.Clustering <- 25
      regressOut_mode <- NULL
      features_to_regressOut <- NULL
      use.sctransform <- TRUE
      vars.to.regress <- c("percent.mt", cell.cycle.features[[regression.mode]])
      cluster.resolution <- 0.5
      
      if (file.exists(file.path(path.to.06.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.rds")) == FALSE){
        DefaultAssay(s.obj) <- "RNA"
        s.obj <- JoinLayers(s.obj)
        s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                             save.RDS.s8 = TRUE,
                                                             path.to.output = path.to.06.output,
                                                             use.sctransform = TRUE,
                                                             num.PCA = num.PCA,
                                                             num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                             num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                             cluster.resolution = cluster.resolution,
                                                             vars.to.regress = vars.to.regress)        
      } else {
        print("integrated done. file at %s", file.path(path.to.06.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.rds"))
      }
    }
  }
}


