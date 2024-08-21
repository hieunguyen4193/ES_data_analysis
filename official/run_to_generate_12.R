#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV5"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "12_process_sub_clusters"
path.to.rmd <- file.path(path.to.main.src, src.dir, "12_downstream_analysis.Rmd")
output_dir <- file.path(path.to.save.html, "12_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

all.filter.cases <- c("BCR_genes",
                      "nCount_and_BCRgenes",
                      "nCount_and_TCRgenes",
                      "nCount",
                      "nCount_and_BCR_TCRgenes",
                      "TCR_genes")

integration.case <- "remove_d4_LPS_SC5"
regression.mode <- "CC_differences"
filter.mode <- "nCount_and_BCR_TCRgenes"

for (subcluster.name in c("b_cells", "lymphoids", "myeloids")){
  output_file <- str_replace(basename(path.to.rmd), ".Rmd", sprintf("%s.%s.%s.%s.html", subcluster.name, integration.case, regression.mode, filter.mode))
  print(sprintf("Working on %s - %s - %s - sub cluster: %s", filter.mode, regression.mode, integration.case, subcluster.name))
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(integration.case = integration.case,
                                    regression.mode = regression.mode,
                                    filter.mode = filter.mode,
                                    subcluster.name = subcluster.name))}      
}

