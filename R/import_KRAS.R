.get_raw_data <- function(genes) {
  raw <- GEOquery::getGEOSuppFiles("GSE137912", baseDir = tempdir(),
                                   filter_regex = "RAW")
  raw_loc <- rownames(raw)
  folder <- stringr::str_remove(raw_loc, "GSE137912_RAW.tar")
  system(paste0("tar -xvf ", raw_loc, " -C ", folder))
  samples <- list.files(folder, full.names = TRUE) %>%
    stringr::str_subset(pattern = "barcodes") %>%
    stringr::str_remove("barcodes.tsv.gz")

  counts_1 <- DropletUtils::read10xCounts(samples[1:4], type = "prefix",
                                          compressed = TRUE)
  keep <- rowData(counts_1)$Symbol %in% genes
  counts_1 <- counts_1[keep, ]
  counts_2 <- DropletUtils::read10xCounts(samples[5:12], type = "prefix",
                                          compressed = TRUE)
  keep <- rowData(counts_2)$Symbol %in% genes
  counts_2 <- counts_2[keep, ]

  keep <- rownames(counts_1) %in% rownames(counts_2)
  counts_1 <- counts_1[keep, ]
  keep <- rownames(counts_2) %in% rownames(counts_1)
  counts_2 <- counts_2[keep, ]
  counts_2 <- counts_2[rownames(counts_1), ]
  colnames(counts_1) <- paste0(counts_1$Sample, counts_1$Barcode)
  colnames(counts_2) <- paste0(counts_2$Sample, counts_2$Barcode)
  counts <- cbind(counts_1, counts_2)
  file.remove(raw_loc)
  rownames(counts) <- rowData(counts)$Symbol
  return(counts(counts))
}

.convert <- function(logcounts) {
  samples <- data.frame(full = colnames(logcounts)) %>%
    dplyr::mutate(root = stringr::word(full, 1, sep = "_"),
                  int = stringr::word(full, 2, sep = "_"),
                  last = stringr::word(full, 3, sep = "_"))
  roots <- samples %>%
    dplyr::select(root, int) %>%
    dplyr::distinct() %>%
    dplyr::group_by(root) %>%
    dplyr::mutate(new_int = rank(int))
  samples_log <- dplyr::full_join(samples, roots) %>%
    dplyr::mutate(new = paste(root, new_int, last, sep = "_")) %>%
    dplyr::select(full, new)
  return(samples_log)
}

.clean_counts <- function(counts) {
  colnames(counts) <- colnames(counts) %>%
    stringr::word(2, sep = "_") %>%
    stringr::str_replace_all("\\.", "_") %>%
    stringr::str_remove("-1$")
  counts <- counts[, !stringr::str_detect(colnames(counts), "H2122_24")]
  samples <- data.frame(full = colnames(counts)) %>%
    dplyr::mutate(root = stringr::word(full, 1, sep = "_"),
                  int = as.numeric(stringr::word(full, 2, sep = "_")),
                  last = stringr::word(full, 3, sep = "_"))
  roots <- samples %>%
    dplyr::select(root, int) %>%
    dplyr::distinct() %>%
    dplyr::group_by(root) %>%
    dplyr::mutate(new_int = rank(int))
  samples_counts <- full_join(samples, roots) %>%
    mutate(new = paste(root, new_int, last, sep = "_"))
  counts <- counts[, samples_counts$full]
  colnames(counts) <- samples_counts$new
  return(counts)
}

#' @title Download and preprocess the raw dataset from Xue et al.
#'
#' @description Download the dataset from GEO, filter, and create a
#' \code{SingleCellExperiment} object
#' @references
#' Jenny Y. Xue, Yulei Zhao, Jordan Aronowitz, Trang T. Mai, Alberto Vides,
#' Besnik Qeriqi, Dongsung Kim, Chuanchuan Li, Elisa de Stanchina, Linas Mazutis,
#'  Davide Risso, and Piro Lito.
#' *Rapid non-uniform adaptation to conformation-specific KRAS(G12C) inhibition.*
#' Nature, 577(7790):421–425,jan  2020. ISSN  14764687. doi: 10.1038/s41586-019-1884-x
#' @return
#' A \code{SingleCellExperiment} object
#' @export
#' @import stringr dplyr readr DropletUtils SingleCellExperiment GEOquery
#' @importFrom openxlsx read.xlsx
import_KRAS <- function() {
  # Get count and logcount matrices ----
  raw <- GEOquery::getGEOSuppFiles("GSE137912", baseDir = tempdir(),
                                   filter_regex = "logcounts")
  raw_loc <- rownames(raw)
  logcounts <- readr::read_csv(raw_loc)
  genes <- logcounts$X1
  logcounts <- logcounts[, -1]
  logcounts <- Matrix::Matrix(as.matrix(logcounts))
  file.remove(raw_loc)
  conversion <- .convert(logcounts)
  counts <- .get_raw_data(genes = genes)
  logcounts <- logcounts[genes %in% rownames(counts), ]
  counts <- .clean_counts(counts)
  counts <- counts[, conversion$new]
  colnames(counts) <- conversion$full
  counts <- counts[rownames(counts) %in% genes, ]
  # Get cell metadata ----
  url <- paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586",
                "-019-1884-x/MediaObjects/41586_2019_1884_MOESM4_ESM.csv")
  celltype <- readr::read_csv(url(url))
  url <- paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586",
                "-019-1884-x/MediaObjects/41586_2019_1884_MOESM6_ESM.csv")
  pst <- readr::read_csv(url(url))
  colD <- full_join(celltype, pst, by = "X1")
  url <- paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586",
                "-019-1884-x/MediaObjects/41586_2019_1884_MOESM13_ESM.xlsx")
  df <- openxlsx::read.xlsx(url, sheet = "Extended Data Fig. 2g, i")
  df <- df[, 1:7]
  colD <- full_join(colD, df)
  rownames(colD) <- colD$X1
  # Get gene metadata ----
  url <- paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586",
                "-019-1884-x/MediaObjects/41586_2019_1884_MOESM5_ESM.csv")
  genemeta <- readr::read_csv(url(url)) %>%
    dplyr::select(-X1, -X22, -Description, -`Gene-level column names`)
  rownames(genemeta) <- genemeta$gene_short_name
  genemeta <- genemeta %>%
    filter(gene_short_name %in% rownames(counts))
  # Format as one object ----
  counts <- counts[genemeta$gene_short_name, ]
  sce <- SingleCellExperiment(assays = list("counts" = counts),
                              rowData = genemeta, colData = DataFrame(colD))
  return(sce)
}
