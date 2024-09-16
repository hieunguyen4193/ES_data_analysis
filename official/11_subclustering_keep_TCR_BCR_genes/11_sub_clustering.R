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

#####-----------------------------------------------------------------------#####
##### SELECTED CONDITIONS
#####-----------------------------------------------------------------------#####
integration.case <- "remove_d4_LPS_SC5"
regression.mode <- "CC_differences"
filter.mode <- "nCount_and_BCR_TCRgenes"

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.10.output <- file.path(path.to.main.output, "10_output", integration.case, regression.mode, filter.mode)
path.to.11.output <- file.path(path.to.main.output, "11_output", integration.case, regression.mode, filter.mode)
dir.create(path.to.11.output, showWarnings = FALSE, recursive = TRUE)
integrate.samples <- integration.config[[integration.case]]

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt", cell.cycle.features[[regression.mode]])
cluster.resolution <- 0.5

s.obj <- readRDS(file.path(path.to.10.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.rds"))

sub.clusters <- list(
  T_cells = c(0, 1, 4, 7, 8, 10, 12, 16, 17),
  B_cells = c(2, 5, 14, 23),
  myeloid = c(3, 6, 9, 11, 13, 15, 16, 18, 19, 20, 21, 22)
)

for (i in names(sub.clusters)){
  if (file.exists(file.path(path.to.11.output, sprintf("subcluster_%s.rds", i))) == FALSE){
    tmp.s.obj <- subset(s.obj, cca.cluster.0.5 %in% sub.clusters[[i]])
    saveRDS(tmp.s.obj, file.path(path.to.11.output, sprintf("subcluster_%s.rds", i)))    
  }
}