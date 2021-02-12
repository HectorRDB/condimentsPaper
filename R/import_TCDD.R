#' @title Download and import the raw dataset from Nault et al.
#' @description Download the dataset from GEO, filter, and create a
#' \code{Seurat object}
#' @references
#' Rance Nault,  Kelly A. Fader,  Sudin Bhattacharya,  and Tim R. Zacharewski.
#' *Single-Nuclei RNA Sequencing Assessment of the Hepatic Effects of*
#' * 2,3,7,8-Tetrachlorodibenzo-p-dioxin.CMGH*,
#' 11(1):147–159, jan 2021.  ISSN 2352345X.
#' doi:  10.1016/j.jcmgh.2020.07.012.
#' @export
#' @import stringr dplyr Seurat GEOquery SingleCellExperiment
#' @importFrom openxlsx read.xlsx
import_TCDD <- function() {
  # Get raw data ----
  raw <- GEOquery::getGEOSuppFiles("GSE148339", baseDir = tempdir(),
                                  filter_regex = "RData")
  raw_loc <- rownames(raw)
  GEOquery::gunzip(raw_loc)
  file <- stringr::str_remove(raw_loc, ".gz")
  cds <- readRDS(file)
  file.remove(file)
  return(cds)
}
