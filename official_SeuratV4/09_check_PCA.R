gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official_SeuratV4"
source(file.path(path.to.project.src, "monocle2_helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"

my_random_seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25

sample1 <- "d4_LPS"
sample2 <- "d4_SPF"

##### change the outdir
# outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
outdir <- "/media/hieunguyen/HD01/outdir/CRC1382"

PROJECT <- "EStange_20240411_SeuratV4"
s.obj.name <- "all_sobj.integrated.rds"
sub.cluster.id <- "subset_231031"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
dir.create(path.to.09.output, showWarnings = FALSE, recursive = TRUE)

path.to.sobj <- file.path(path.to.04.output, s.obj.name, "preprocessed_subcluster_obj", sprintf("preprocessed_subclusters_%s.rds", sub.cluster.id))
s.obj <- readRDS(path.to.sobj)
object.name <- str_replace(basename(path.to.sobj), ".rds", "")
s.obj <- subset(s.obj, name %in% c(sample1, sample2))

