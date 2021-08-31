#' @title Download and preprocess the raw dataset from .
#'
#' @description Download the dataset from GEO, filter, and create a
#' \code{SingleCellExperiment} object
#' @references
#' Arun C. Habermann et al.
#' *Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis*
#' Nature, 577(7790):421–425,jan  2020. ISSN  14764687. doi: 10.1038/s41586-019-1884-x
#' @return
#' A \code{SingleCellExperiment} object
#' @export
#' @import Seurat SingleCellExperiment
import_fibrosis <- function() {
  # Get Seurat Object
  url <- ""
  cds <- readRDS(file = url(url))
  cds <- cds[, cds$population == "Epithelial"]
  cds <- cds[, cds$celltype %in% c("AT1", "AT2", "SCGB3A2+", "Transitional AT2")]
  sce <- as.SingleCellExperiment(cds)
  filt <- apply(counts(sce), 1, function(x){
    sum(x >= 2) >= 30
  })
  assays(sce) <- assays(sce)[1]
  altExp(sce) <- NULL
  sce <- sce[filt, ]
  return(sce)
}
