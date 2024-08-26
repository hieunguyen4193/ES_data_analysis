#####-----------------------------------------------------------------------#####
##### RUN TO GENERATE 01 HTML REPORTS
#####-----------------------------------------------------------------------#####

gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
# PROJECT <- "EStange_20240411_SeuratV5"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

src.dir <- "01_QC_hashtag_antibodies_data"

all.samples <- c("adult_GF",
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
                 "SC11")

#####----------------------------------------------------------------------#####
##### generate html for 01_count_cells_with_hashtags.Rmd
#####----------------------------------------------------------------------#####
for (sample.id in all.samples){
  if (sample.id == "d7_GF"){
    path.to.rmd <- file.path(path.to.main.src, src.dir, "01_count_cells_with_hashtags.d7_GF.Rmd")  
  } else if (sample.id == "SC11"){
    print(path.to.rmd)
    path.to.rmd <- file.path(path.to.main.src, src.dir, "01_count_cells_with_hashtags.SC11.Rmd")
  } else {
    path.to.rmd <- file.path(path.to.main.src, src.dir, "01_count_cells_with_hashtags.Rmd")
  }
  print(sample.id)
  print(path.to.rmd)
  output_file <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", sample.id))
  output_dir <- file.path(path.to.save.html, "01_output")
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(sample.id = sample.id, 
                                    PROJECT = PROJECT, 
                                    outdir = outdir))    
  }
}

#####----------------------------------------------------------------------#####
##### generate html for 01_comparing_different_quantile.Rmd
#####----------------------------------------------------------------------#####
for (sample.id in all.samples){
  if (sample.id != "d7_GF"){
    path.to.rmd <- file.path(path.to.main.src, src.dir, "01_comparing_different_quantile.Rmd")  
  } else {
    path.to.rmd <- file.path(path.to.main.src, src.dir, "01_comparing_different_quantile.d7_GF.Rmd")
  }
  output_file <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", sample.id))
  output_dir <- file.path(path.to.save.html, "01_output")
  if (file.exists(file.path(output_dir, output_file)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = output_file,
                      output_dir = output_dir,
                      params = list(sample.id = sample.id, 
                                    PROJECT = PROJECT, 
                                    outdir = outdir))    
  }
}
