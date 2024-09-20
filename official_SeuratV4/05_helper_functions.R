run_pseudobulk_dge <- function(input.s.obj, sample1, sample2){
  s.obj <- hash()
  s.obj[[sample1]] <- subset(input.s.obj, name == sample1)
  s.obj[[sample2]] <- subset(input.s.obj, name == sample2)
  
  DefaultAssay(s.obj[[sample1]]) <- "RNA"
  DefaultAssay(s.obj[[sample2]]) <- "RNA"
  
  hashtag.cells <- hash()
  hashtag.cells[[sample1]] <- hash()
  hashtag.cells[[sample2]] <- hash()
  
  all.ht <- hash()
  
  all.ht[[sample1]] <- setdiff(unique(s.obj[[sample1]]$HTO_classification), c("Negative"))
  all.ht[[sample2]] <- setdiff(unique(s.obj[[sample2]]$HTO_classification), c("Negative"))
  
  for (sample.id in c(sample1, sample2)){
    true.exists.ht <- c()
    for (ht in all.ht[[sample.id]]){
      ht.cells <- colnames(subset(s.obj[[sample.id]], HTO_classification == ht))
      if (length(ht.cells) > 1){
        hashtag.cells[[sample.id]][[ht]] <- ht.cells      
        true.exists.ht <- c(true.exists.ht, ht)
      }
    }
    all.ht[[sample.id]] <- true.exists.ht
  }
  
  count.matrix <- hash()
  
  pseudobulk.matrix <- data.frame(Gene = rownames(s.obj[[sample1]]))
  for (sample.id in c(sample1, sample2)){
    count.matrix[[sample.id]] <- hash()
    tmp.count.matrix <- GetAssayData(object = s.obj[[sample.id]], slot = "counts", assay = "RNA")
    for (ht in all.ht[[sample.id]]){
      tmp <- tmp.count.matrix[, hashtag.cells[[sample.id]][[ht]]]
      rowsum.tmp <- rowSums(tmp)
      tmpdf <- data.frame(data = rowsum.tmp[rowsum.tmp >= 10]) %>% rownames_to_column("Gene")
      colnames(tmpdf) <- c("Gene", sprintf("%s_%s", sample.id, ht))
      count.matrix[[sample.id]][[ht]] <- tmpdf
      pseudobulk.matrix <- merge(pseudobulk.matrix, tmpdf, by.x = "Gene", by.y = "Gene")
    }
  }
  
  pseudobulk.matrix <- pseudobulk.matrix %>% column_to_rownames("Gene")
  pseudo_metadata <- data.frame(sample.id = colnames(pseudobulk.matrix))
  pseudo_metadata$group <- unlist(lapply(pseudo_metadata$sample.id, function(x){
    str_split(x, "_Hashtag")[[1]][[1]]
  }))
  
  pseudo_metadata$group <- factor(pseudo_metadata$group, levels = c(sample1, sample2))
  
  dds <- DESeqDataSetFromMatrix(pseudobulk.matrix, 
                                colData = pseudo_metadata, 
                                design = ~ group)
  
  dds <- DESeq(dds)
  contrast <- c("group", levels(pseudo_metadata$group)[2], levels(pseudo_metadata$group)[1])
  
  # resultsNames(ddsrow
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05)
  
  resdf <- res %>%
    data.frame() %>% rownames_to_column("Gene")
  
  padj_cutoff <- 0.05
  logFC_cutoff <- 1
  
  normalized_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("Gene")
  resdf <- merge(resdf, normalized_counts, by.x = "Gene", by.y = "Gene") %>%
    mutate(abslog2FoldChange = abs(log2FoldChange)) %>%
    mutate(sig = ifelse(padj <= padj_cutoff & abslog2FoldChange >= logFC_cutoff, "sig.diff", "nonsig.diff")) 
  
  resdf.sigOnly <- subset(resdf, resdf$padj <= padj_cutoff & resdf$abslog2FoldChange >= logFC_cutoff) %>% arrange(desc(log2FoldChange))
  if (nrow(resdf.sigOnly) == 0){
    heatmap.full.input.scaled <- data.frame()
    heatmap.top50.input <- data.frame()
  } else {
    heatmap.input <- rbind(head(resdf.sigOnly, 10), tail(resdf.sigOnly, 10))[, c("Gene", pseudo_metadata$sample.id)]
    
    heatmap.top50.input <- rbind(head(resdf.sigOnly, 50), tail(resdf.sigOnly, 50))[, c("Gene", pseudo_metadata$sample.id)]
    
    heatmap.full.input <-resdf.sigOnly[, c("Gene", pseudo_metadata$sample.id)] %>% column_to_rownames("Gene")
    heatmap.full.input.scaled <- (heatmap.full.input - rowSums(heatmap.full.input))/rowSds(as.matrix(heatmap.full.input)) 
    
    heatmap.full.input.scaled <- heatmap.full.input.scaled %>% rownames_to_column("Gene") 
    heatmap.full.input.scaled <- merge(heatmap.full.input.scaled, subset(resdf.sigOnly, select = c(Gene, padj, log2FoldChange, abslog2FoldChange)),
                                       by.x = "Gene", by.y = "Gene") %>% arrange(desc(log2FoldChange))
    
  }
  output <- list(
    dds = dds,
    resdf = resdf,
    resdf.sigOnly = resdf.sigOnly,
    heatmap.top50.input = heatmap.top50.input,
    heatmap.full.input.scaled = heatmap.full.input.scaled,
    pseudo_metadata = pseudo_metadata
  )
  return(output)
}