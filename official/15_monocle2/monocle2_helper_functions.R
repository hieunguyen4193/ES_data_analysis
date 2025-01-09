#####
##### RUN MONOCLE2 FROM S.OBJ
#####
run_monocle2_from_presave_obj <- function(monocle.obj, path.to.save.monocle.obj){
  library(monocle)
  print("Estimate size factors ...")
  monocle.obj <- estimateSizeFactors(monocle.obj)
  print("Estimate dispersions ...")
  monocle.obj <- estimateDispersions(monocle.obj)
  
  print("Detect expressed genes ...")
  monocle.obj <- detectGenes(monocle.obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle.obj),
                                      num_cells_expressed >= 100))
  
  fData(monocle.obj)$use_for_ordering <-
    fData(monocle.obj)$num_cells_expressed > 0.05 * ncol(monocle.obj)
  
  ordering.genes <- subset(fData(monocle.obj), fData(monocle.obj)$use_for_ordering == TRUE)$gene_short_name
  ordering.genes <- intersect(ordering.genes, expressed_genes)
  monocle.obj <- monocle.obj[ordering.genes,]
  
  set.seed(my_random_seed)
  print("Running dimensional reduction ...")
  monocle.obj <- reduceDimension(monocle.obj,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 3,
                                 reduction_method = 'tSNE',
                                 verbose = T, 
                                 random_seed = my_random_seed)
  set.seed(my_random_seed)
  
  print("running cluster cells ...")
  monocle.obj <- clusterCells(monocle.obj, 
                              verbose = TRUE, 
                              random_seed = my_random_seed)
  
  print("running clustering DEG genes ...")
  clustering_DEG_genes <-
    differentialGeneTest(monocle.obj[ordering.genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 35, verbose = TRUE)
  
  HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  print("Set ordering filter...")
  monocle.obj <-
    setOrderingFilter(monocle.obj,
                      ordering_genes = HSMM_ordering_genes)
  set.seed(my_random_seed)
  
  print("Running dimensional reduction with DDRTree algorithm ...")
  monocle.obj <-
    reduceDimension(monocle.obj, method = 'DDRTree', random_seed = my_random_seed)
  set.seed(my_random_seed)
  
  print("Running orderCells...")
  monocle.obj <- orderCells(monocle.obj)
  
  print("plot cell tracjectory...")
  p <- plot_cell_trajectory(monocle.obj)
  
  print("saving monocle2 results to rds object")
  saveRDS(monocle.obj, file.path(path.to.save.monocle.obj, sprintf("monocle_obj.rds")))
  print("finished!")
  return(monocle.obj)
}
