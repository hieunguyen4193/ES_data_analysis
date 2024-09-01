#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 04 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4_reduced_RNAcontam_0"
s.obj.name <- "all_sobj.integrated.rds"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

path.to.rmd <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4/05_DEA_between_d4_LPS_vs_d4_SPF.Rmd"

output_dir <- file.path(path.to.save.html, "05_output")

for (sub.cluster.id in c("Myeloid_Basophils", "B_cells", "T_cells")){
  output_file <- sprintf("05_DEA_between_d4_LPS_vs_d4_SPF_%s.html", sub.cluster.id)
  
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(sub.cluster.id = sub.cluster.id, 
                                    PROJECT = PROJECT, 
                                    outdir = outdir,
                                    s.obj.name = s.obj.name))    
  }
}
