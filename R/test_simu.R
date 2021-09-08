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
  da_results <- da_results %>% dplyr::slice_min(SpatialFDR, n = 1, with_ties = FALSE)
  return(da_results)
}


.running_slingshot <- function(sce, shape = 2) {
  if (shape == 5) {
    ends <- c('sD', 'sF', 'sG', 'sH', 'sI')
    clusters <- str_remove(sce$from, "mid")
    sds <- slingshot(reducedDim(sce), clusters, start.clus = "sA",
                     end.clus = ends)
    
  } else {
    if (shape == 2) {
      ends <- c("sC", 'sD')
    } else {
      ends <- c('sE', "sF", "sG")
    }
    sds <- slingshot(reducedDim(sce), sce$from, start.clus = "sA",
                     end.clus = ends, approx_points = 100, stretch = 1)
  }
  return(sds)
}


.running_monocle <- function(sce, clusters, start = 1,
                             params = list(orthogonal_proj_tip = TRUE),
                             keep_cds = FALSE) {
  colnames(sce) <- paste0("Cell-", seq_len(ncol(sce)))
  names(clusters) <- colnames(sce)
  fd <- data.frame(gene_short_name = rownames(sce))
  rownames(fd) <- rownames(sce)
  pd <- data.frame(cellid = colnames(sce))
  rownames(pd) <- colnames(sce)
  cds <- new_cell_data_set(counts(sce), cell_metadata = pd, gene_metadata = fd)
  cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds)
  cds@int_colData$reducedDims$UMAP <- reducedDim(sce[,colnames(cds)])
  cds <- cluster_cells(cds)
  partitions <- cds@clusters$UMAP$clusters <- clusters[colnames(cds)]
  partitions[colnames(cds)] <- 1
  cds@clusters$UMAP$partitions <- partitions %>% 
    as.factor() %>%
    droplevels()
  cds <- learn_graph(cds, learn_graph_control = params,
                     close_loop = FALSE)
  cds <- order_cells(cds, root_cells = colnames(cds)[clusters == start])
  plot_cells(cds)
  mst <- principal_graph(cds)$UMAP
  roots <- cds@principal_graph_aux$UMAP$root_pr_nodes
  # Get the other endpoints
  endpoints <- names(which(igraph::degree(mst) == 1))
  root <- endpoints[endpoints %in% roots]
  endpoints <- endpoints[!endpoints %in% roots]
  # See https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/#subset-branch
  cellWeights <- lapply(endpoints, function(endpoint) {
    # We find the path between the endpoint and the root
    path <- choose_graph_segments(cds, ending_pr_nodes = endpoint, 
                                  starting_pr_node = root, return_list = TRUE)
    # We find the cells that map along that path
    df <- data.frame(weights = as.numeric(colnames(cds) %in% path$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .) %>%
    as.matrix()
  rownames(cellWeights) <- colnames(cds)
  pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
                       nrow = ncol(cds), byrow = FALSE)
  colnames(pseudotime) <- colnames(cellWeights)
  res <- list("pseudotime" = pseudotime,
              "cellWeights" = cellWeights,
              "conditions" = sce[, colnames(cds)]$condition)
  if(keep_cds) {
    res$cds <- cds
  }
  return(res)
}


.condiments_analysis <- function(sce, shape = 2) {
  sds <- running_slingshot(sce, shape)
  # clusters
  clusters <- sce$from %>% 
    str_remove("mid") %>% 
    as.factor() %>% 
    as.numeric() %>% 
    as.factor()
  names(clusters) <- colnames(sce)
  # monocle
  cds <- running_monocle(sce, clusters)
  vals <- bind_rows(
    "sling_prog" = progressionTest(sds, sce$condition, thresh = .01, 
                                   lineages = TRUE, global = FALSE),
    "mon_prog" = progressionTest(pseudotime = cds$pseudotime, 
                                 cellWeights = cds$cellWeights,
                                 conditions = cds$condition, thresh = .01,
                                 lineages = TRUE, global = FALSE),
    "sling_diff" = fateSelectionTest(sds, sce$condition, method = "Classifier",
                                     args_classifier = list(method = "rf"), 
                                     thresh = .01, pairwise = TRUE, global = FALSE),
    .id = "test_type"
  ) %>%
    group_by(test_type) %>%
    mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    dplyr::slice_min(FDR, n = 1, with_ties = FALSE) %>%
    arrange(test_type)
  
  vals$nLineages <- c(ncol(cds$cellWeights),
                      nLineages(sds),
                      nLineages(sds)
  )
  return(vals)
}

.condiments_analysis_per_cond <- function(sce, shape = 2) {
  # Sling
  sds <- running_slingshot(sce, shape)
  sdss <- slingshot_conditions(sds, sce$condition)
  n <- nLineages(sdss[[1]])
  sds <- merge_sds(sdss[[1]], sdss[[2]], condition_id = names(sdss),
                   mapping = matrix(1:n, nrow = n, ncol = 2))
  vals <- bind_rows(
    "sling_prog" = progressionTest(sds, sce$condition, thresh = .05, 
                                   lineages = TRUE, global = FALSE),
    "sling_diff" = fateSelectionTest(sds, sce$condition,
                                     method = "Classifier",
                                     args_classifier = list(method = "rf"),
                                     thresh = .01, 
                                     pairwise = TRUE, global = FALSE),
    .id = "test_type"
  ) %>%
    group_by(test_type) %>%
    mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    dplyr::slice_min(FDR, n = 1, with_ties = FALSE)
  vals$nLineages <- c(mean(nLineages(sdss[[1]]), nLineages(sdss[[2]])),
                      mean(nLineages(sdss[[1]]), nLineages(sdss[[2]]))
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
    colMedians(reducedDim(sce)[labels == label, ])  %>% matrix(nrow = 1) %>%
      return()
  }) %>% t()
  da_regions$FDR <- spatialFDR(intensity, da_regions$pval.wilcoxon)
  da_regions <- da_regions %>%
    dplyr::rename("statistic" = "DA.score",
                  "p.value" = "pval.wilcoxon") %>%
    mutate(test_type = paste0("region", rownames(.))) %>%
    select(test_type, FDR, p.value, statistic)
  return(da_regions %>% dplyr::slice_min(FDR, n = 1, with_ties = FALSE))
}

#' @title Analyse dataset with two conditions
#' @param sce The dataset, a SingleCellExperiment object
#' @param shape Number of lineages (2, 3 or 5)
#' @export
#' @import miloR DAseq slingshot dplyr SingleCellExperiment monocle3
#' @importFrom cydar spatialFDR
#' @importFrom matrixStats colMedians
anayze_all <- function(sce, shape = 2) {
  res <- tryCatch({suppressMessages(
    suppressWarnings(
      bind_rows(
        "condiments" = .condiments_analysis(sce, shape = shape) %>%
          select(-lineage, -pair),
        "DAseq" = .DAseq_analysis(sce),
        "milo" = .milo_analysis(sce) %>%
          select(PValue, FDR, SpatialFDR, `F`) %>%
          dplyr::rename("p.value" = PValue,
                        "statistic" = `F`),
        .id = "method"
      ) %>% as.data.frame() %>%
        return()
    )
  )}, error = function(cond) {
    return(sce)
  })
  return(res)
}

#' @title Analyse dataset with three conditions
#' @param sce The dataset, a SingleCellExperiment object
#' @param shape Number of lineages (2 or 3)
#' @export
#' @import miloR slingshot dplyr SingleCellExperiment monocle3
anayze_multiple_conditions <- function(sce) {
  res <- tryCatch({suppressMessages(
    suppressWarnings(
      bind_rows(
        "condiments" = .condiments_analysis(sce, shape = shape) %>%
          select(-lineage, -pair),
        "milo" = .milo_analysis(sce) %>%
          select(PValue, FDR, SpatialFDR, `F`) %>%
          dplyr::rename("p.value" = PValue,
                        "statistic" = `F`),
        .id = "method"
      )
    )
  )}, error = function(cond) {
    return(sce)
  })
  return(res)
}

#' @title Analyse dataset with common and separate trajectories per condition
#' @param sce The dataset, a SingleCellExperiment object
#' @param shape Number of lineages (2, 3 or 5.)
#' @export
#' @import slingshot dplyr SingleCellExperiment monocle3
anayze_all_per_cond <- function(sce, shape = 2) {
  res <- tryCatch({suppressMessages(
    suppressWarnings(
      bind_rows(
        "normal" = .condiments_analysis(sce, shape = shape) %>%
          select(-lineage, -pair),
        "failure" = .condiments_analysis_per_cond(sce, shape = shape) %>%
          select(-lineage, -pair),
        .id = "method"
      )
    )
  )}, error = function(cond) {
    print(cond)
    return(sce)
  })
  return(res)
}