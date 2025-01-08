gc()
rm(list = ls())

library(stringdist)
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

# path.to.save.topClone.fasta <- file.path(path.to.main.output, "topClone_FASTA_20250107")
path.to.save.topClone.fasta <- file.path(path.to.main.output, "topClone_FASTA_20250107_full")
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

##### update 07.01.2025
##### subset: keep only Hashtag1
# s.obj <- subset(s.obj, HTO_classification == "Hashtag1-TotalSeqC")

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% subset(select = -c(CTaa))
all.clonedf <- subset(all.clonedf, all.clonedf$barcode %in% colnames(s.obj)) 

# merge clone information with the metadata of the seurat object. 
all.clonedf <- merge(all.clonedf, meta.data, by.x = "barcode", by.y = "barcode")

for (chain in c("TRA", "TRB")){
  input.clonedf <- subset(all.clonedf, all.clonedf$chain == chain)
  dir.create(file.path(path.to.save.topClone.fasta, chain), showWarnings = FALSE, recursive = TRUE)
  
  trab.countdf <- table(input.clonedf$cdr3) %>% as.data.frame()
  colnames(trab.countdf) <- c("CTaa", "count")
  trab.countdf <- trab.countdf %>% arrange(desc(count)) %>% head(300)
  trab.countdf$cloneID <- to_vec(
    for (item in seq(1, 300)) sprintf("ClonalType%s", item)
  )
  
  seqs <- c()
  min.dist <- c()
  max.dist <- c()
  for (x in trab.countdf$CTaa){
    tmpdf <- subset(input.clonedf, input.clonedf$cdr3 == x) %>%
      rowwise() %>%
      mutate(seq = paste0(substr(fwr3_nt, nchar(fwr3_nt) - 20 + 1, nchar(fwr3_nt)),
                          cdr3_nt, 
                          substr(fwr4_nt, 1, 20)))
    count.tmpdf <- table(tmpdf$seq) %>% data.frame() %>% arrange(desc(Freq))
    selected.seq <- subset(count.tmpdf, count.tmpdf$Freq == max(count.tmpdf$Freq))$Var1[[1]]
    seqs <- c(seqs, as.character(selected.seq))
    
    count.tmpdf <- count.tmpdf %>% rowwise() %>%
      mutate(SELECTED = ifelse(Var1 == selected.seq, "YES", "NO")) %>%
      mutate(dist.to.selected.seq = stringdist( Var1, selected.seq, method="lv"))
    if (nrow(count.tmpdf) > 1){
      min.dist <- c(min.dist, min(subset(count.tmpdf, count.tmpdf$SELECTED == "NO")$dist.to.selected.seq))
      max.dist <- c(max.dist, min(subset(count.tmpdf, count.tmpdf$SELECTED == "NO")$dist.to.selected.seq))
    } else {
      min.dist <- c(min.dist, "unique sequence")
      max.dist <- c(max.dist, "unique sequence")
    }
    
    
    writexl::write_xlsx(count.tmpdf, file.path(path.to.save.topClone.fasta, chain, sprintf("%s_count_and_edit_distance.xlsx", x)))
  }
  trab.countdf$seq <- seqs
  trab.countdf$min.dist <- min.dist
  trab.countdf$max.dist <- max.dist
  writexl::write_xlsx(trab.countdf, file.path(path.to.save.topClone.fasta, sprintf("%s.top300.xlsx", chain)))
  
  path.to.output.fasta <- file.path(path.to.save.topClone.fasta, 
                                    sprintf("%s.fasta", chain))
  ##### save to FASTA files
  sink(path.to.output.fasta)
  for (i in seq(1, nrow(trab.countdf))){
    clonenal.type <- trab.countdf[i, ]$cloneID
    output.seq <- trab.countdf[i, ]$CDR3seq
    output.info <- sprintf(">%s", clonenal.type)            
    cat(output.info)
    cat("\n")
    cat(output.seq)
    cat("\n")
  }
  sink()
}


