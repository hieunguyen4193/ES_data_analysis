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

input.ht <- "Hashtag1-TotalSeqC"
# input.ht <- "Hashtag2-TotalSeqC"

# path.to.save.topClone.fasta <- file.path(path.to.main.output, "topClone_FASTA_20250107")
path.to.save.topClone.fasta <- file.path(path.to.main.output, "topclone_FASTA_20250114", input.ht)
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
s.obj <- subset(s.obj, HTO_classification == input.ht)

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% subset(select = -c(CTaa))
all.clonedf <- subset(all.clonedf, all.clonedf$barcode %in% colnames(s.obj)) 

# merge clone information with the metadata of the seurat object. 
all.clonedf <- merge(all.clonedf, meta.data, by.x = "barcode", by.y = "barcode")

add.seq <- "CAAERNSNNRIFF_CAWSRTGEDTQYF"
add.seqs <- list(TRA = "CAAERNSNNRIFF",
                 TRB = "CAWSRTGEDTQYF")

for (input.chain in c("TRA", "TRB")){
  input.clonedf <- subset(all.clonedf, all.clonedf$chain == input.chain)
  dir.create(file.path(path.to.save.topClone.fasta, input.chain), showWarnings = FALSE, recursive = TRUE)
  
  trab.countdf <- table(input.clonedf$cdr3) %>% as.data.frame()
  colnames(trab.countdf) <- c("CTaa", "count")
  trab.countdf <- trab.countdf %>% arrange(desc(count)) %>% head(300)
  trab.countdf$cloneID <- to_vec(
    for (item in seq(1, 300)) sprintf("ClonalType%s", item)
  )
  if (input.ht == "Hashtag1-TotalSeqC"){
    trab.countdf <- rbind(
      trab.countdf, 
      data.frame(CTaa = add.seqs[[input.chain]], 
                 count = nrow(subset(all.clonedf, all.clonedf$cdr3 == add.seqs[[input.chain]])),
                 cloneID = "Add.ClonalType")
    )    
  }
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
    
    
    writexl::write_xlsx(count.tmpdf, file.path(path.to.save.topClone.fasta, input.chain, sprintf("%s_count_and_edit_distance.xlsx", x)))
  }
  trab.countdf$seq <- seqs
  trab.countdf$min.dist <- min.dist
  trab.countdf$max.dist <- max.dist
  writexl::write_xlsx(trab.countdf, file.path(path.to.save.topClone.fasta, sprintf("%s.top300.xlsx", input.chain)))
  
  path.to.output.fasta <- file.path(path.to.save.topClone.fasta, 
                                    sprintf("%s.fasta", input.chain))
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

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("barcode") %>% 
  subset(select = c(barcode, CTaa)) %>%
  rowwise() %>%
  subset(is.na(CTaa) == FALSE) %>%
  mutate(input.barcode = barcode) %>%
  subset(grepl("NA_", CTaa) == FALSE) %>%
  subset(grepl("_NA", CTaa) == FALSE) %>%
  mutate(TRA = str_split(CTaa, "_")[[1]][[1]]) %>%
  mutate(TRB = str_split(CTaa, "_")[[1]][[2]]) %>%
  mutate(TRA_V_gene = paste(subset(all.clonedf, all.clonedf$chain == "TRA" & all.clonedf$barcode == input.barcode)$v_gene, collapse = ",")) %>%
  mutate(TRA_J_gene = paste(subset(all.clonedf, all.clonedf$chain == "TRA" & all.clonedf$barcode == input.barcode)$j_gene, collapse = ",")) %>%
  mutate(TRB_V_gene = paste(subset(all.clonedf, all.clonedf$chain == "TRB" & all.clonedf$barcode == input.barcode)$v_gene, collapse = ",")) %>%
  mutate(TRB_J_gene = paste(subset(all.clonedf, all.clonedf$chain == "TRB" & all.clonedf$barcode == input.barcode)$j_gene, collapse = ","))

meta.data <- meta.data %>% rowwise() %>%
  mutate(VJ_CDR3_TRA = sprintf("%s_%s_%s", TRA_V_gene, TRA_J_gene, TRA)) %>%
  mutate(VJ_CDR3_TRB = sprintf("%s_%s_%s", TRB_V_gene, TRB_J_gene, TRB))

countdf <- data.frame( num.unique.CTaa = c(length(unique(meta.data$CTaa))),
                       num.unique.TRA = c(length(unique(meta.data$TRA))),
                       num.unique.TRB = c(length(unique(meta.data$TRB))),
                       num.unique.TRA_VJ = c(length(unique(meta.data$VJ_CDR3_TRA))),
                       num.unique.TRB_VJ = c(length(unique(meta.data$VJ_CDR3_TRB))),
                       dupl.TRA = nrow(meta.data[duplicated(meta.data$TRA), ]),
                       dupl.TRB = nrow(meta.data[duplicated(meta.data$TRB), ]))

writexl::write_xlsx(countdf, file.path(path.to.save.topClone.fasta, "count_unique_clone.xlsx"))
