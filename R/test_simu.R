.milo_analysis <- function(sce) {
  logcounts(sce) <- log1p(counts(sce))
  milo <- Milo(sce)
  milo <- buildGraph(milo, reduced.dim = "dimred", d = 3, BPPARAM = SerialParam(),
                     k = 20)
  milo <- makeNhoods(milo, refined = TRUE, reduced_dims = "dimred", d = 3)
  # We simulate replicate in each condition, similar to
  # https://github.com/MarioniLab/milo_analysis_2020/blob/a84308a01675c4b4daae0fddddbd1ada70679f9b/notebooks/SFig2_batch_effect_simulation.Rmd#L52
  milo$Sample <- sample(1:5, ncol(sce), replace = TRUE)
  for (i in seq_along(unique(sce$condition))) {
    cond <- unique(sce$condition)[i]
    milo$Sample[milo$condition == cond] <-
      milo$Sample[sce$condition == cond] + 5 * (i - 1)
  }
  milo$Sample <- as.character(milo$Sample)
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)),
                            sample = "Sample")
  milo <- calcNhoodDistance(x = milo, d = 3, reduced.dim = "dimred")
  da_results <- testNhoods(milo, design = ~ condition,
                           design.df = colData(milo) %>%
                             as.data.frame() %>%
                             select(condition, Sample) %>%
                             distinct())
  da_results <- da_results %>% dplyr::slice_min(SpatialFDR, 1, with_ties = FALSE)
  return(da_results)
}

.running_slingshot <- function(sce, shape = 2) {
  if (shape == 2) {
    ends <- c("sC", 'sD')
  } else {
    ends <- c('sE', "sF", "sG")
  }

  sds <- slingshot(reducedDim(sce), sce$from, start.clus = "sA",
                   end.clus = ends, approx_points = 100, stretch = 1)
  return(sds)
}

.condiments_analysis <- function(sce, shape = 2) {
  sds <- .running_slingshot(sce, shape)
  vals <- bind_rows(
    "prog" = progressionTest(sds, sce$condition, thresh = .01),
    "diff" = differentiationTest(sds, sce$condition, method = "Classifier",
                                 args_classifier = list(method = "rf"), thresh = .01),
    .id = "test_type"
  )
  return(vals)
}

.DAseq_analysis <- function(sce) {
  colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  rownames(reducedDim(sce)) <- colnames(sce)
  sce$Sample <- sample(1:10, ncol(sce), replace = TRUE)
  sce$Sample[sce$condition == "WT"] <- sce$Sample[sce$condition == "WT"] + 10
  # sce$Sample <- sce$condition
  sce$Sample <- as.character(sce$Sample)
  info <- colData(sce) %>%
    as.data.frame() %>%
    select(Sample, condition) %>%
    distinct()
  da_cells <- getDAcells(
    X = reducedDim(sce, "dimred"),
    cell.labels = sce$Sample,
    labels.1 = info$Sample[info$condition == "WT"],
    labels.2 = info$Sample[info$condition == "KO"],
    k.vector = seq(50, 500, 50),
    do.plot = FALSE
  )
  if(length(da_cells$da.up) < 2 | length(da_cells$da.down) < 2) {
    return(data.frame("test_type" = "region1", "FDR" = 1,
                      "p.value" = 1, "statistic" = 0))
  }
  da_regions <- tryCatch({
    getDAregion(
      X = reducedDim(sce, "dimred"),
      da.cells = da_cells,
      cell.labels = sce$Sample,
      labels.1 = info$Sample[info$condition == "WT"],
      labels.2 = info$Sample[info$condition == "KO"],
      do.plot = FALSE
    )}, error = function(cond) {
      return(list("DA.stat" = data.frame()))
    })
  labels <- da_regions$da.region.label
  da_regions <- da_regions$DA.stat %>% as.data.frame()
  if (nrow(da_regions) == 0) {
    return(
      data.frame("test_type" = "region1", "FDR" = 1,
                 "p.value" = 1, "statistic" = 0))
  }
  if (nrow(da_regions) == 1) {
    da_regions$FDR <- da_regions$pval.wilcoxon
    return(da_regions)
  }
  intensity <- sapply(rownames(da_regions) %>% as.numeric, function(label) {
    matrixStats::colMedians(reducedDim(sce)[labels == label, ]) %>%
      matrix(nrow = 1) %>%
      return()
  }) %>% t()
  da_regions$FDR <- cydar::spatialFDR(intensity, da_regions$pval.wilcoxon)
  da_regions <- da_regions %>%
    dplyr::rename("statistic" = "DA.score",
                  "p.value" = "pval.wilcoxon") %>%
    mutate(test_type = paste0("region", rownames(.))) %>%
    select(test_type, FDR, p.value, statistic)
  return(da_regions %>% dplyr::slice_min(FDR, 1, with_ties = FALSE))
}


#' @title Analyse dataset with two conditions
#' @param sce The dataset, a SingleCellExperiment object
#' @param shape Number of lineages (2 or 3)
#' @export
#' @import miloR DAseq slingshot dplyr SingleCellExperiment
#' @importFrom cydar spatialFDR
#' @importFrom matrixStats colMedians
anayze_all <- function(sce, shape = 2) {
  res <- tryCatch({suppressMessages(
    suppressWarnings(
      bind_rows(
        "condiments" = .condiments_analysis(sce, shape = 2) %>%
          select(-lineage, -pair),
        "DAseq" = .DAseq_analysis(sce),
        "milo" = .milo_analysis(sce) %>%
          dplyr::select(PValue, FDR, SpatialFDR, `F`) %>%
          dplyr::rename("p.value" = PValue,
                        "statistic" = `F`),
        .id = "method"
      )
    )
  )}, error = function(cond) {
    print(cond)
    return(sce)
  })
  return(res)
}

#' @title Analyse dataset with three conditions
#' @param sce The dataset, a SingleCellExperiment object
#' @param shape Number of lineages (2 or 3)
#' @export
#' @import miloR slingshot dplyr SingleCellExperiment
anayze_multiple_conditions <- function(sce) {
  res <- tryCatch({suppressMessages(
    suppressWarnings(
      bind_rows(
        "condiments" = .condiments_analysis(sce, shape = 2) %>%
          select(-lineage, -pair),
        "milo" = .milo_analysis(sce) %>%
          dplyr::select(PValue, FDR, SpatialFDR, `F`) %>%
          dplyr::rename("p.value" = PValue,
                        "statistic" = `F`),
        .id = "method"
      )
    )
  )}, error = function(cond) {
    print(cond)
    return(sce)
  })
  return(res)
}

