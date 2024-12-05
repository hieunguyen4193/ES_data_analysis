
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  
  if (require("data.table", quietly = TRUE)) {
    data.table::fwrite(
      data.table::data.table(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv")
    )
  } else {
    write.table(
      data.frame(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv"),
      sep = ",", row.names = FALSE
    )
  }
  
  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
  
  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )
  
  if (file.info(bundle)$size / 1e6 > 500) {
    warning("The output file is more than 500 MB and will need to be subset further.")
  }
}

writeCounts <- function(out_dir, data, barcodes = colnames(data), gene.id = rownames(data), gene.symbol = rownames(data), feature.type = "Gene Expression", subset = 1:length(barcodes)) {
  require("R.utils")
  require("Matrix")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (require("data.table")) {
    data.table::fwrite(
      data.table::data.table(barcode = barcodes[subset]),
      file.path(out_dir, "barcodes.tsv.gz"),
      col.names = FALSE
    )
    
    data.table::fwrite(
      data.table::data.table(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      col.names = FALSE
    )
  } else {
    write.table(
      data.frame(barcode = barcodes[subset]),
      gzfile(file.path(out_dir, "barcodes.tsv.gz")),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
    
    write.table(
      data.frame(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
  }
  
  Matrix::writeMM(data[, subset], file.path(out_dir, "matrix.mtx"))
  R.utils::gzip(file.path(out_dir, "matrix.mtx"), remove = TRUE)
}