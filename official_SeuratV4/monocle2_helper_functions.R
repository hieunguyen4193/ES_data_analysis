
#####
##### RUN MONOCLE2 FROM S.OBJ
#####
run_monocle2 <- function(s.obj, path.to.save.monocle.obj){
  library(monocle)
  data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
  
  pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
  
  fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fd)
  
  monocle.obj <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  monocle.obj <- estimateSizeFactors(monocle.obj)
  monocle.obj <- estimateDispersions(monocle.obj)
  
  monocle.obj <- detectGenes(monocle.obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle.obj),
                                      num_cells_expressed >= 100))
  
  fData(monocle.obj)$use_for_ordering <-
    fData(monocle.obj)$num_cells_expressed > 0.05 * ncol(monocle.obj)
  
  ordering.genes <- subset(fData(monocle.obj), fData(monocle.obj)$use_for_ordering == TRUE)$gene_short_name
  ordering.genes <- intersect(ordering.genes, expressed_genes)
  monocle.obj <- monocle.obj[ordering.genes,]
  
  set.seed(my_random_seed)
  monocle.obj <- reduceDimension(monocle.obj,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 3,
                                 reduction_method = 'tSNE',
                                 verbose = T, 
                                 random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- clusterCells(monocle.obj, verbose = F, random_seed = my_random_seed)
  
  clustering_DEG_genes <-
    differentialGeneTest(monocle.obj[ordering.genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 20)
  
  HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  monocle.obj <-
    setOrderingFilter(monocle.obj,
                      ordering_genes = HSMM_ordering_genes)
  set.seed(my_random_seed)
  monocle.obj <-
    reduceDimension(monocle.obj, method = 'DDRTree', random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- orderCells(monocle.obj)
  
  p <- plot_cell_trajectory(monocle.obj)
  
  saveRDS(monocle.obj, file.path(path.to.save.monocle.obj, sprintf("monocle_obj.rds")))
  return(monocle.obj)
}

#####----------------------------------------------------------------------#####
##### extract monocle information from monocle object, helper function
#####----------------------------------------------------------------------#####
extract_monocle_info <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    stop(paste0("For now tradeSeq only support Monocle with DDRTree",
                "reduction. If you want to use another type",
                "please use another format for tradeSeq inputs."))
  }
  # Get the reduced dimension of DDRT
  rd <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
  
  # Get the various lineages info for weights and pseudotime
  y_to_cells <- cds@auxOrderingData[["DDRTree"]]
  y_to_cells <- y_to_cells$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  root <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root <- y_to_cells$Y[y_to_cells$cells == root]
  mst <- monocle::minSpanningTree(cds)
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[endpoints != paste0("Y_", root)]
  cellWeights <- lapply(endpoints, function(endpoint) {
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .)
  pseudotime <- sapply(cellWeights, function(w) cds$Pseudotime)
  rownames(cellWeights) <- rownames(pseudotime) <- colnames(cds)
  # Get the lineages representation
  edges_rd <- t(monocle::reducedDimK(cds)) %>% as.data.frame()
  rd_lineages <- lapply(endpoints, function(endpoint){
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    path <- paste("Y", path, sep = "_")
    return(edges_rd[path, ])
  })
  return(list("pseudotime" = pseudotime,
              "cellWeights" = as.matrix(cellWeights)))
}


#####----------------------------------------------------------------------#####
##### Get earliest principal node 
#####----------------------------------------------------------------------#####
get_earliest_principal_node <- function(cds, cluster.id){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == cluster.id)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}