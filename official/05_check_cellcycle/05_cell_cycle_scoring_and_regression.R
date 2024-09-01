gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/EStange/official"
source(file.path(path.to.project.src, "project.config.R"))

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "EStange_20240411_reduced_RNAcontam_0"
  
stage_lst <- hash()
stage_lst[["adult_GF"]] <- c(adult_GF = "adult_GF")
stage_lst[["adult_SPF"]] <- c(adult_SPF = "adult_SPF")
stage_lst[["d4_GF"]] <- c(d4_GF = "d4_GF")
stage_lst[["d4_LPS"]] <- c(d4_LPS = "d4_LPS")
stage_lst[["d4_SPF"]] <- c(d4_SPF = "d4_SPF")
stage_lst[["d7_GF"]] <- c(d7_GF = "d7_GF")
stage_lst[["d7_SPF"]] <- c(d7_SPF = "d7_SPF")
stage_lst[["d10_SPF"]] <- c(d10_SPF = "d10_SPF")
stage_lst[["d15_SPF"]] <- c(d15_SPF = "d15_SPF")
stage_lst[["d20_SPF"]] <- c(d20_SPF = "d20_SPF")
stage_lst[["SC5"]] <- c(SC5 = "SC5")
stage_lst[["SC11"]] <- c(SC11 = "SC11")
stage_lst[["SC12"]] <- c(SC12 = "SC12")

for (sample.id in names(stage_lst)){
  for (cellcycle.regression.mode in c("S_G2M_G1_scores", "CC_differences")){
    print (sprintf("Working on sample %s, mode: %s", sample.id, cellcycle.regression.mode))
    path.to.main.input <- file.path(outdir, PROJECT)
    path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
    path.to.01.output <- file.path(path.to.main.output, "01_output")
    path.to.02.output <- file.path(path.to.main.output, "02_output", sample.id)
    path.to.05.output <- file.path(path.to.main.output, "05_output", sample.id, cellcycle.regression.mode)
    dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
    
    s.obj <- readRDS(file.path(path.to.02.output, sprintf("GEX_sample_%s_seurat_object.rds", sample.id)))
    
    all.genes <- rownames(x = s.obj)
    s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
    s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
    g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
    g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
    s.obj <- CellCycleScoring(s.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    s.obj$G1.Score = 1 - s.obj$S.Score - s.obj$G2M.Score
    s.obj$CC.Difference <- s.obj$S.Score - s.obj$G2M.Score
    
    s.obj <- RunPCA(object = s.obj, features = c(s.genes, g2m.genes), reduction.name = "CELLCYCLE_PCA_before")
    cellcycle.pca.before <- DimPlot(object = s.obj, reduction = "CELLCYCLE_PCA_before", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "Phase")
    
    umap.cellcycle.score <- FeaturePlot(object = s.obj, reduction = "RNA_UMAP", features = c("S.Score", "G2M.Score", "G1.Score"))
    umap.cellcycle.phase <- DimPlot(object = s.obj, reduction = "RNA_UMAP", group.by = "Phase", label = TRUE, label.box = TRUE)
    
    input.cluster <- 1
    cellcycle.marker.plot <- RidgePlot(subset(s.obj, seurat_clusters == input.cluster  ), features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
    
    if (cellcycle.regression.mode == "S_G2M_G1_scores"){
      s.obj <- SCTransform(object = s.obj, assay = "SCT", new.assay.name = "SCT", vars.to.regress = c("S.Score", "G2M.Score", "G1.Score"))  
    } else if (cellcycle.regression.mode == "CC_differences"){
      s.obj <- SCTransform(object = s.obj, assay = "SCT", new.assay.name = "SCT", vars.to.regress = c("CC.Difference"))  
    }
    
    s.obj <- RunPCA(object = s.obj, features = c(s.genes, g2m.genes), reduction.name = "CELLCYCLE_PCA_after")
    cellcycle.pca.after <- DimPlot(object = s.obj, reduction = "CELLCYCLE_PCA_after", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "Phase")
    
    compare.pca.cellcycle <- cellcycle.pca.before + cellcycle.pca.after
    saveRDS(object = s.obj, file.path(path.to.05.output, sprintf("%s.cellcycle_reg.rds", sample.id)))
    ggsave(plot = compare.pca.cellcycle, filename = "compare_PCA_before_vs_after_cellcycle_regression.svg", path = path.to.05.output, device = "svg", width = 14, height = 10, dpi = 300)
  }
}



