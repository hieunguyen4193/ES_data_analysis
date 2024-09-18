run_integration_after_subsetting <- function(input.s.obj){
  my_random_seed <- 42
  
  num.dim.integration <- 25 
  num.PCA <- 25
  num.dim.cluster <- 25
  num.PC.used.in.Clustering <- 25
  num.PC.used.in.UMAP <- 25
  cluster.resolution  <- 0.5
  
  input.s.obj <- DietSeurat(input.s.obj)
  
  chosen.assay <- "RNA"
  DefaultAssay(input.s.obj) <- chosen.assay
  
  ##### Re run UMAP and integration
  input.s.obj <- NormalizeData(input.s.obj) # ---> use Log Normalized
  input.s.obj <- FindVariableFeatures(input.s.obj, selection.method = "vst")
  input.s.obj <- ScaleData(input.s.obj, features = rownames(input.s.obj))
  
  input.s.obj <- RunPCA(input.s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  
  input.s.obj <- RunUMAP(input.s.obj, reduction = sprintf("%s_PCA", chosen.assay), 
                         dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_UMAP", chosen.assay),
                         seed.use = my_random_seed, umap.method = "uwot")
  
  # clustering 
  input.s.obj <- FindNeighbors(input.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  input.s.obj <- FindClusters(input.s.obj, resolution = cluster.resolution, random.seed = 0)
  
  #### Integration
  data.list <- SplitObject(input.s.obj, split.by = "name")
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
  
  k.filter <- 200
  
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                    k.filter = k.filter) ## THIS IS CCA DIMENSIONS
  
  input.s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter) ## THIS IS PCA DIMENSION
  
  ## keep the order of integration obj
  input.s.obj_inte <- input.s.obj_inte[, colnames(input.s.obj)]
  
  input.s.obj[['integrated']] <- input.s.obj_inte[['integrated']]
  
  input.s.obj@commands <- c(input.s.obj@commands, input.s.obj_inte@commands)
  
  input.s.obj@tools <- c(input.s.obj@tools, input.s.obj_inte@tools)
  
  DefaultAssay(input.s.obj) <- "integrated"
  input.s.obj <- ScaleData(input.s.obj, verbose = FALSE, features = row.names(input.s.obj))
  
  input.s.obj <- RunPCA(input.s.obj, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA")
  
  input.s.obj <- RunUMAP(input.s.obj, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP")
  
  # clustering 
  input.s.obj <- FindNeighbors(input.s.obj, reduction = "INTE_PCA", dims = 1:num.dim.cluster)
  
  input.s.obj <- FindClusters(input.s.obj, resolution = cluster.resolution)
  return(input.s.obj)
}
