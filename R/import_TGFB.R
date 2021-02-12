#' @title Download and import the raw dataset from McFaline et al..
#'
#' @description Download the dataset from GEO, filter, and create a
#' \code{SingleCellExperiment object}
#' @references
#' Jos\'{e}  L.  McFaline-Figueroa,  Andrew  J.  Hill,  Xiaojie  Qiu,
#' Dana  Jackson,  Jay  Shendure,  and  Cole Trapnell.
#' *A pooled single-cell genetic screen identifies regulatory checkpoints in *
#' *the continuum of the epithelial-to-mesenchymal transition.*
#' Nature Genetics, 51(9):1389–1398, sep 2019.  ISSN 15461718.
#' doi:  10.1038/s41588-019-0489-5.
#' @export
#' @import SingleCellExperiment slingshot GEOquery
import_TGFB <- function(){
  raw <- GEOquery::getGEOSuppFiles("GSE114687", baseDir = tempdir(),
                                   filter_regex = "GSE114687_pseudospace_cds")
  raw_loc <- rownames(raw)
  GEOquery::gunzip(raw_loc)
  file <- stringr::str_remove(raw_loc, ".gz")
  cds <- readRDS(file)
  file.remove(file)

  # Extract useful info from the cellDataSet object
  counts <- cds@assayData$exprs
  phenoData <- pData(cds@phenoData)
  rd <- SimpleList(
    tSNEorig = cbind(cds@phenoData@data$TSNE.1, cds@phenoData@data$TSNE.2)
  )
  filt <- apply(counts, 1, function(x){
    sum(x >= 2) >= 15
  })
  counts <- counts[filt, ]
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = phenoData,
    reducedDims = rd)
  return(sce)
}
