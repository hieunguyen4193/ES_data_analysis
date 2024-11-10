gc()
rm(list = ls())

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

path.to.main.output <- file.path(outdir, PROJECT)
path.to.save.html <- file.path(path.to.main.output, "html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"

integration.case <- "remove_d4_LPS_SC5"
regression.mode <- "CC_differences"
filter.mode <- "nCount_and_BCR_TCRgenes"
subcluster.name <- "T_cells"
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.12.output <- file.path(path.to.main.output, "12_output_remove_BCR_TCR", integration.case, regression.mode, filter.mode, subcluster.name)

path.to.save.topClone.fasta <- file.path(path.to.main.output, "topClone_FASTA")
dir.create(path.to.save.topClone.fasta, showWarnings = FALSE, recursive = TRUE)

path.to.VDJ.input <- "/media/hieunguyen/HD01/storage/EStange/EStange_20240411/VDJ"
path.to.VDJ.output <- file.path(outdir, PROJECT, "VDJ_output")
# all.clone.files <- Sys.glob(file.path(path.to.VDJ.output, "annotated_contigs_clonaltype_*.csv"))
all.clone.files <- Sys.glob(file.path(path.to.VDJ.input, "*", "filtered_contig_annotations.csv"))

all.clonedf <- data.frame()
for (input.file in all.clone.files){
  sampleid <- dirname(input.file) %>% basename()
  tmpdf <- read.csv(input.file) %>%
    rowwise() %>%
    mutate(barcode = sprintf("%s_%s", sampleid, barcode))
  all.clonedf <- rbind(all.clonedf, tmpdf)
}

s.obj <- readRDS(file.path(path.to.12.output, "s8_output", "EStange_20240411_reduced_RNAcontam_0.output.s8.rds"))
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% subset(select = -c(CTaa))
all.clonedf <- subset(all.clonedf, all.clonedf$barcode %in% colnames(s.obj)) 

# merge clone information with the metadata of the seurat object. 
all.clonedf <- merge(all.clonedf, meta.data, by.x = "barcode", by.y = "barcode")

#####----------------------------------------------------------------------#####
##### Split the clone information to IGH and IGL
#####----------------------------------------------------------------------#####

tra.clonedf <- subset(all.clonedf, all.clonedf$chain == "TRA")
trb.clonedf <- subset(all.clonedf, all.clonedf$chain == "TRB")

tra.countdf <- table(tra.clonedf$cdr3) %>% as.data.frame()
colnames(tra.countdf) <- c("CTaa", "count")
tra.countdf <- tra.countdf %>% arrange(desc(count)) %>% head(300)
tra.countdf$cloneID <- to_vec(
  for (item in seq(1, 300)) sprintf("ClonalType%s", item)
)

tra.countdf$CDR3seq <- unlist(lapply(
  tra.countdf$CTaa, function(x){
    tmpdf <- subset(tra.clonedf, tra.clonedf$cdr3 == x)
    tmp.count.seq <- table(tmpdf$cdr3_nt)
    max.tmp.count.seq <- tmp.count.seq[tmp.count.seq == max(tmp.count.seq)] %>% names()
    return(max.tmp.count.seq[[1]] )
  }
))

trb.countdf <- table(trb.clonedf$cdr3) %>% as.data.frame()
colnames(trb.countdf) <- c("CTaa", "count")
trb.countdf <- trb.countdf %>% arrange(desc(count)) %>% head(300)
trb.countdf$cloneID <- to_vec(
  for (item in seq(1, 300)) sprintf("ClonalType%s", item)
)

trb.countdf$CDR3seq <- unlist(lapply(
  trb.countdf$CTaa, function(x){
    tmpdf <- subset(trb.clonedf, trb.clonedf$cdr3 == x)
    tmp.count.seq <- table(tmpdf$cdr3_nt)
    max.tmp.count.seq <- tmp.count.seq[tmp.count.seq == max(tmp.count.seq)] %>% names()
    return(max.tmp.count.seq[[1]] )
  }
))

path.to.output.fasta <- file.path(path.to.save.topClone.fasta, "TRA.fasta")
##### save to FASTA files
sink(path.to.output.fasta)
for (i in seq(1, nrow(tra.countdf))){
  clonenal.type <- tra.countdf[i, ]$cloneID
  output.seq <- tra.countdf[i, ]$CDR3seq
  output.info <- sprintf(">%s", clonenal.type)            
  cat(output.info)
  cat("\n")
  cat(output.seq)
  cat("\n")
}
sink()

path.to.output.fasta <- file.path(path.to.save.topClone.fasta, "TRB.fasta")
##### save to FASTA files
sink(path.to.output.fasta)
for (i in seq(1, nrow(trb.countdf))){
  clonenal.type <- trb.countdf[i, ]$cloneID
  output.seq <- trb.countdf[i, ]$CDR3seq
  output.info <- sprintf(">%s", clonenal.type)            
  cat(output.info)
  cat("\n")
  cat(output.seq)
  cat("\n")
}
sink()

writexl::write_xlsx(tra.countdf, file.path(path.to.save.topClone.fasta, "TRA.top300.xlsx"))
writexl::write_xlsx(trb.countdf, file.path(path.to.save.topClone.fasta, "TRB.top300.xlsx"))