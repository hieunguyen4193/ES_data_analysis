#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 02 HTML REPORTS
#####-----------------------------------------------------------------------#####
gc()
rm(list = ls())

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "EStange_20240411_SeuratV4"
config.version <- "reduced_RNAcontam_0"

PROJECT <- sprintf("%s_%s", PROJECT, config.version)

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

all.samples <- c("adult_GF",
                 "d4_GF",   
                 "adult_SPF",
                 "d4_LPS",
                 "d10_SPF",
                 "d4_SPF",
                 "d15_SPF",
                 "d20_SPF",
                 "d7_SPF",
                 "SC12",
                 "SC11",
                 "d7_GF", 
                 "SC5")

##### install clusterProfiler
if ("clusterProfiler" %in% installed.packages() == FALSE){
  remove.packages("clusterProfiler")
  remove.packages("DOSE")
  remove.packages("GOSemSim")
  path.to.install.dir <- "/media/hieunguyen/HD01/storage/offline_pkgs/clusterProfiler"
  install.packages(file.path(path.to.install.dir, "HDO.db_0.99.1.tar.gz"), type = "source", repos = NULL)
  install.packages(file.path(path.to.install.dir, "yulab.utils_0.1.4.tar.gz"), type = "source", repos = NULL)
  install.packages(file.path(path.to.install.dir, "GOSemSim_2.28.1.tar.gz"), type = "source", repos = NULL)
  install.packages(file.path(path.to.install.dir, "DOSE_3.28.2.tar.gz"), type = "source", repos = NULL) 
  install.packages(file.path(path.to.install.dir, "gson_0.1.0.tar.gz"), type = "source", repos = NULL)
  install.packages(file.path(path.to.install.dir, "clusterProfiler_4.10.1.tar.gz"), type = "source", repos = NULL) 
}

for (sample.id in all.samples){
  if (sample.id == "SC11"){
    path.to.rmd <- file.path(path.to.main.src, "02_preliminary_analysis", "02_GEX_preliminary_data_analysis.SC11.Rmd")
    output_file <- str_replace(basename(path.to.rmd), ".SC11.Rmd", sprintf(".%s.html", sample.id))
  } else if (sample.id == "SC5"){
    path.to.rmd <- file.path(path.to.main.src, "02_preliminary_analysis", "02_GEX_preliminary_data_analysis.SC5.Rmd")
    output_file <- str_replace(basename(path.to.rmd), ".SC5.Rmd", sprintf(".%s.html", sample.id))
  } else {
    path.to.rmd <- file.path(path.to.main.src, "02_preliminary_analysis", "02_GEX_preliminary_data_analysis.Rmd")
    output_file <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", sample.id))
  }
  print(sample.id)
  print(path.to.rmd)
  
  output_dir <- file.path(path.to.save.html, "02_output")
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(sample.id = sample.id, 
                                    PROJECT = PROJECT, 
                                    outdir = outdir))    
  }
}