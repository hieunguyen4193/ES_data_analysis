gc()
rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggrepel)

s.obj <- readRDS("/home/hieunguyen/CRC1382/outdir/tmp/preprocessed_subclusters_subset_231031_JW.rds")

p.original <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)

# convert.cluster.to.name <- list(0 = "celltype1",
#                                 2 = "celltype2")

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(new.cluster = ifelse(seurat_clusters %in% c(0, 1, 2), 0, seurat_clusters)) %>%
  # mutate(new.cluster = convert.cluster.to.name[[seurat_clusters]]) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$new.cluster, col.name = "new.cluster")
  
p.new <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "new.cluster")

s.obj.d4 <- subset(s.obj, name %in% c("d4_SPF", "d4_LPS"))

p.d4 <- DimPlot(object = s.obj.d4, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "new.cluster")

p.d4.name <- DimPlot(object = s.obj.d4, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, split.by = "name")

DimPlot(object = s.obj.d4, reduction = "INTE_PCA", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name")
DimPlot(object = s.obj.d4, reduction = "RNA_PCA", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name")

#####----------------------------------------------------------------------#####
##### Subset sample d4 only
#####----------------------------------------------------------------------#####
list.of.genes <- c("Elane", "Mpo", "Chil3","Ly6c2", "Mgst1","Fn1", "Sell", "F13a1", "Cd177","Ace", "Gngt2", "Adgre4", "Itgal", "Cx3cr1", "Tmem176b", "Slamf7", "H2-Aa")
s.obj.d4 <- RunPCA(object = s.obj.d4, assay = "RNA", npcs = 2, features = list.of.genes, reduction.name = "test.PCA")
DimPlot(object = s.obj.d4, reduction = "test.PCA", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name")

#####----------------------------------------------------------------------#####
##### Subset cluster 11 only
#####----------------------------------------------------------------------#####
s.obj.11 <- subset(s.obj.d4, seurat_clusters == 11)
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
cluster.resolution <- 0.5
pca_reduction_name <- "PCA.11"
umap_reduction_name <- "UMAP.11"
my_random_seed <- 42

#####----------------------------------------------------------------------#####
##### Re-run PCA and UMAP after subsetting cluster 11
#####----------------------------------------------------------------------#####
s.obj.11 <- RunPCA(s.obj.11, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
s.obj.11 <- RunUMAP(s.obj.11, reduction = pca_reduction_name, 
                    dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                    seed.use = my_random_seed, umap.method = "uwot")

#####----------------------------------------------------------------------#####
##### Run re-clustering after subsetting the data at cluster 11
#####----------------------------------------------------------------------#####
s.obj.11 <- FindNeighbors(s.obj.11, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
s.obj.11 <- FindClusters(s.obj.11, resolution = cluster.resolution, random.seed = 0)
cluster.markers11 <- FindAllMarkers(object = s.obj.11, assay = "RNA", test.use = "wilcox")
cluster.markers11 <- subset(cluster.markers11, cluster.markers11$p_val_adj <= 0.05 & cluster.markers11$avg_log2FC > 0)

#####----------------------------------------------------------------------#####
##### UMAP gene expression plot, pick top N genes.
#####----------------------------------------------------------------------#####
plot.list <- list()
for (input.cluster in unique(cluster.markers11$cluster)){
  tmp.cluster.markers <- subset(cluster.markers11, cluster.markers11$cluster == input.cluster) %>% arrange(desc(avg_log2FC)) 
  DefaultAssay(s.obj.11) <- "RNA"
  plot.list[[sprintf("cluster_%s", input.cluster)]] <-  FeaturePlot(object = s.obj.11, reduction = "UMAP.11", 
                                                                    features = head(tmp.cluster.markers, 12)$gene, 
                                                                    label = TRUE, 
                                                                    repel = TRUE, 
                                                                    ncol = 4)  
  # plot.list$cluster_... <<<
}

p.cluster11.name <- DimPlot(object = s.obj.11, reduction = "UMAP.11", label = TRUE, group.by = "name")

#####----------------------------------------------------------------------#####
##### Diff. test between conditions, volcano plot
#####----------------------------------------------------------------------#####
# input.cluster <- 0
# subset(s.obj.11, seurat_clusters == input.cluster)

diff.genes <- FindMarkers(object = s.obj.11, group.by = "name", ident.1 = "d4_LPS", ident.2 = "d4_SPF")
cutoff.adjp <- 0.05

diff.genes <- diff.genes %>% 
  rownames_to_column("Gene") %>%
  rowwise() %>%
  mutate(sig = ifelse(p_val_adj <= 0.05 & abs(avg_log2FC) > 1, "sig.", "not sig.")) %>%
  mutate(show.gene.name = ifelse(sig == "sig.", Gene, NA))



volcano.plot <- ggplot(data=diff.genes, 
                       aes(x=avg_log2FC, y=-log10(p_val_adj), col = sig, label=Gene)) + 
  geom_point() + geom_label_repel(label = diff.genes$show.gene.name, size = 8) + 
  scale_color_manual(values=c("#c0d2f0", "#f28095")) +
  theme_minimal() +
  geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
  geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
  #geom_text() +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12), axis.text=element_text(size=12))

p.violin.plot1 <- VlnPlot(object = s.obj.11, features = c("Sept5"), split.by = "name", split.plot = TRUE)

input.genes <- c("Sept5", "Itgb2l", "4930438A08Rik", "Stfa2l1")
count.matrix <- GetAssayData(object = s.obj.11, assay = "RNA", slot = "data")[input.genes, ] %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode")
count.matrix.pivot <- count.matrix %>% pivot_longer(!barcode, names_to = "Gene", values_to = "exprs") 
count.matrix.pivot <- merge(as.data.frame(count.matrix.pivot), s.obj.11@meta.data %>% rownames_to_column("barcode") ,
                            by.x = "barcode", by.y = "barcode")

# count.matrix.pivot %>% ggplot(aes(x = Gene, y = exprs, fill = name)) + geom_boxplot() 
count.matrix.pivot %>% ggplot(aes(x = Gene, y = exprs, fill = name)) + geom_violin()
