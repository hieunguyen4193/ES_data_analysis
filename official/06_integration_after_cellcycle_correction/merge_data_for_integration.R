gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.R"))

options(future.globals.maxSize = 10000 * 1024^2)

library(argparse)

# Create a parser object
parser <- ArgumentParser()

# Define arguments
parser$add_argument("--integration_case", required = TRUE, help = "Choose the integration case")
parser$add_argument("--regression_mode", required = TRUE, help = "Choose the regression mode")

# Parse the arguments
args <- parser$parse_args()

# Access the arguments
regression.mode <- args$regression_mode
integration.case <- args$integration_case


# for i in all.samples remove_d4_LPS_SC5 remove_d4_LPS remove_d4_LPS_SC5_SC11;do \
# for j in CC_differences S_G2M_G1_scores;do Rscript merge_data_for_integration.R \ 
# --integration_case $i --regression_mode $j;done;done

#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####
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

# for (regression.mode in c("CC_differences", "S_G2M_G1_scores")){
#   for (integration.case in names(integration.config)){
    print(sprintf("working on integration %s with regression mode %s", integration.case, regression.mode))
    integrate.samples <- integration.config[[integration.case]]
    ##### merging the data before integration
    if (file.exists(file.path(path.to.merge.data, sprintf("finished_saving_data_case_%s_%s.csv", integration.case, regression.mode))) == FALSE){
      print("Merge data does not exists, generate new merged data...")
      data.list <- list()
      for (i in seq(length(integrate.samples))){
        sample.id <- integrate.samples[[i]]
        path.to.05.output <- file.path(path.to.main.output, "05_output", sample.id, regression.mode)
        print(sprintf("reading in sample %s", sample.id))
        data.list[[i]] <- readRDS(file.path(path.to.05.output, sprintf("%s.cellcycle_reg.rds", sample.id)))
      }
      s.obj <- merge(data.list[[1]], data.list[2: length(integrate.samples)])
      saveRDS(s.obj, file.path(path.to.merge.data, sprintf("raw_merge_dataset_%s_%s.rds", integration.case, regression.mode)))  
      write.csv(data.frame(status = c("finished saving data")), 
                file.path(path.to.merge.data, sprintf("finished_saving_data_case_%s_%s.csv", integration.case, regression.mode)))
      } 
#   }
# }


