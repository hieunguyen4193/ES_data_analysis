gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5/processes_src"

source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering_SeuratV5.selectedGenes.R"))

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

integrate.samples <- integration.config[[integration.case]]

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt", cell.cycle.features[[regression.mode]])
cluster.resolution <- 0.5
subcluster.name <- "T_cells"
path.to.12.output <- file.path(path.to.main.output, 
                               "12_output_remove_BCR_TCR", 
                               integration.case, 
                               regression.mode, 
                               filter.mode, 
                               subcluster.name)

Tcell.subcluster.version <- "v1"
path.to.17.output <- file.path(path.to.main.output, 
                               "17_output_remove_BCR_TCR", 
                               integration.case, 
                               regression.mode, 
                               filter.mode, 
                               subcluster.name, 
                               Tcell.subcluster.version)
dir.create(path.to.17.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(
  path.to.12.output, 
  "s8_output", 
  "EStange_20240411_reduced_RNAcontam_0.output.s8.rds"
))

# Eva's message on 29.10.2024
all.Tcell.sub.clusters <- list(
  v1 = c(0,1,3,4,5,7,9,10,11,12,13)
)
selected.CD4.clusters <- all.Tcell.sub.clusters[[Tcell.subcluster.version]] 

subset.s.obj <- subset(s.obj, cca.cluster.0.5 %in% selected.CD4.clusters)
all.genes <- row.names(s.obj)
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
nCount <- s.obj@meta.data$nCount_RNA
tmp.metadata <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  mutate(complexity = log10(nFeature_RNA) / log10(nCount_RNA) )
all.qt <- quantile(nCount) %>% as.list()
iqr <- all.qt$`75%` - all.qt$`25%`
upper.cut <- all.qt$`75%` + 2.5 * iqr
selected.cells <- subset(tmp.metadata, tmp.metadata$nCount_RNA <= upper.cut & tmp.metadata$complexity >= 0.8)

if (filter.mode == "TCR_genes"){
  remove.genes <- TCRgenes.to.exclude
} else if (filter.mode == "nCount"){
  s.obj <- subset(s.obj, cells = selected.cells$barcode)
  remove.genes <- NULL
} else if (filter.mode == "BCR_genes"){
  remove.genes <- BCRgenes.to.exclude      
} else if (filter.mode == "nCount_and_TCRgenes"){
  s.obj <- subset(s.obj, cells = selected.cells$barcode)
  remove.genes <- TCRgenes.to.exclude
} else if (filter.mode == "nCount_and_BCRgenes"){
  s.obj <- subset(s.obj, cells = selected.cells$barcode)
  remove.genes <- BCRgenes.to.exclude
} else if (filter.mode == "BCR_TCR_genes"){
  remove.genes <- c(TCRgenes.to.exclude, BCRgenes.to.exclude)
} else if (filter.mode == "nCount_and_BCR_TCRgenes"){
  s.obj <- subset(s.obj, cells = selected.cells$barcode)
  remove.genes <- c(TCRgenes.to.exclude, BCRgenes.to.exclude)
}

s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = subset.s.obj, 
                                                     save.RDS.s8 = TRUE,
                                                     path.to.output = path.to.17.output,
                                                     use.sctransform = TRUE,
                                                     num.PCA = num.PCA,
                                                     num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                     num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                     cluster.resolution = cluster.resolution,
                                                     vars.to.regress = vars.to.regress,
                                                     remove.genes = remove.genes)        

write.csv(data.frame(status = c("finished saving...")), 
          file.path(path.to.17.output, "s8_output", "EStange_20240411_SeuratV5.output.s8.csv"))

  