.calculate_dimred <- function (model, num_landmarks = 1000, dimred_simulations = TRUE, 
                              dimred_gold = TRUE) {
  is_tf <- NULL
  prev_num_cores <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = model$num_cores)
  sim_counts <- model$simulations$counts
  sim_ix <- seq_len(nrow(sim_counts))
  gs_counts <- model$gold_standard$counts
  gs_ix <- length(sim_ix) + seq_len(nrow(gs_counts))
  tf_info <- model$feature_info %>% filter(is_tf)
  sim_wcounts <- sim_counts[, tf_info$mol_premrna, drop = FALSE]
  sim_xcounts <- sim_counts[, tf_info$mol_mrna, drop = FALSE]
  sim_ycounts <- sim_counts[, tf_info$mol_protein, drop = FALSE]
  gs_wcounts <- gs_counts[, tf_info$mol_premrna, drop = FALSE]
  gs_xcounts <- gs_counts[, tf_info$mol_mrna, drop = FALSE]
  gs_ycounts <- gs_counts[, tf_info$mol_protein, drop = FALSE]
  counts <- rbind(sim_xcounts, gs_xcounts)
  counts@x <- log2(counts@x + 1)
  max_cols <- apply(counts, 2, quantile, 0.99)
  max_cols[max_cols == 0] <- 1
  counts <- sweep(counts, 2, max_cols, "/") %>% Matrix::Matrix(sparse = TRUE)
  landmark_ix <- c(gs_ix, sim_ix)
  if (length(landmark_ix) > num_landmarks) {
    landmark_ix <- sample(landmark_ix, num_landmarks)
  }
  dist_2lm <- as.matrix(
    dynutils::calculate_distance(x = counts[landmark_ix, 
                                            , drop = FALSE],
                                 y = counts, method = model$distance_metric))
  attr(dist_2lm, "landmark_ix") <- landmark_ix
  dimred <- lmds::cmdscale_landmarks(dist_2lm, ndim = 4)
  dimnames(dimred) <- list(rownames(counts), 
                           paste0("comp_", seq_len(ncol(dimred))))
  model$simulations$dimred <- dimred[sim_ix, , drop = FALSE]
  model$gold_standard$dimred <- dimred[gs_ix, , drop = FALSE]
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = prev_num_cores)
  return(model)
}

# Inspired from
# https://github.com/dynverse/dyngen/tree/devel#vignettes
.convert_to_sce <- function(model) {
  counts <-  model$simulations$counts
  counts <- counts[, str_detect(colnames(counts), "mol_mrna")]
  counts <- counts %>% as.matrix() %>% t()
  rownames(counts) <- str_remove(rownames(counts), "mol_mrna_")
  meta <-  model$simulations$meta %>%
    mutate(condition = case_when(
      str_detect(from, "WT") ~ "WT",
      str_detect(from, "UP") ~ "UP",
      TRUE ~ "KO"),
      from = word(from, 2, sep = "_"),
      to = word(to, 2, sep = "_"))
  dimred <-  model$simulations$dimred
  feature <- model$feature_info
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = meta,
                              rowData = feature)
  reducedDim(sce, "dimred") <- dimred
  return(sce)
}

.generate_model_common <- function(backbone, nSim) {
  model_common <-
    initialise_model(
      backbone = backbone,
      num_cells = 100,
      num_tfs = nrow(backbone$module_info),
      num_targets = 250,
      num_hks = 250,
      simulation_params = simulation_default(
        census_interval = 10,
        ssa_algorithm = ssa_etl(tau = 300 / 3600),
        experiment_params = simulation_type_wild_type(num_simulations = nSim)
      )
    ) %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard()
  return(model_common)
}

.generate_model_comb <- function(model_common, module, multiplier, nSim) {
  model_wt <- model_common %>%
    generate_cells()
  b3_genes <- model_common$feature_info %>% filter(module_id == module) %>% pull(feature_id)
  model_ko <- model_common
  model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = nSim,
    timepoint = 0,
    genes = b3_genes,
    num_genes = length(b3_genes),
    multiplier = multiplier
  )
  model_ko <- model_ko %>%
    generate_cells()
  model_comb <-
    combine_models(list(WT = model_wt, KO = model_ko)) %>%
    generate_experiment()
  return(model_comb)
}

.generate_model_three <- function(model_common, multiplier, nSim) {
  model_wt <- model_common %>%
    generate_cells()
  b3_genes <- model_common$feature_info %>% filter(module_id == "B3") %>% pull(feature_id)
  model_ko <- model_common
  model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = nSim,
    timepoint = 0,
    genes = b3_genes,
    num_genes = length(b3_genes),
    multiplier = multiplier
  )
  model_ko <- model_ko %>%
    generate_cells()
  model_up <- model_common
  model_up$simulation_params$experiment_params <- simulation_type_knockdown(
    num_simulations = nSim,
    timepoint = 0,
    genes = b3_genes,
    num_genes = length(b3_genes),
    multiplier = round(1 / multiplier, 2)
  )
  model_up <- model_up %>%
    generate_cells()
  model_comb <-
    combine_models(list(WT = model_wt, KO = model_ko, UP = model_up)) %>%
    generate_experiment()
  return(model_comb)
}

#' @title Create a bifurcating trajectory with two conditions
#' @param multiplier Effect on the knockout
#' @param nSim number of simulations for each condition
#' @export
#' @import SingleCellExperiment dyngen stringr dplyr
create_bifurcating_simu <- function(multiplier = 1, nSim = 200) {
  backbone <- backbone_bifurcating()
  model_common <- .generate_model_common(backbone, nSim)
  model_comb <- .generate_model_comb(model_common, "B3", multiplier, nSim)
  sce <- .convert_to_sce(model_comb)
  return(sce)
}

#' @title Create a trajectory with two bifurcations and two conditions
#' @param multiplier Effect on the knockout
#' @param nSim number of simulations for each condition
#' @export
#' @import SingleCellExperiment dyngen stringr dplyr
create_consecutive_bifurcating_simu <- function(multiplier = 1, nSim = 200) {
  backbone <- backbone_consecutive_bifurcating()
  model_common <- .generate_model_common(backbone, nSim)
  model_comb <- .generate_model_comb(model_common, "D2", multiplier, nSim)
  sce <- .convert_to_sce(model_comb)
  return(sce)
}

#' @title Create a bifurcating trajectory with three conditions
#' @param multiplier Effect on the knockout
#' @param nSim number of simulations for each condition
#' @export
#' @import SingleCellExperiment dyngen stringr dplyr
create_bifurcating_three_conditions_simu <- function(multiplier = 1, nSim = 200) {
  backbone <- backbone_bifurcating()
  model_common <- .generate_model_common(backbone, nSim)
  model_comb <- .generate_model_three(model_common, multiplier, nSim)
  sce <- .convert_to_sce(model_comb)
  return(sce)
}

#' @title Create a bifurcating trajectory with five lineages and three conditions
#' @param multiplier Effect on the knockout
#' @param nSim number of simulations for each condition
#' @export
#' @import SingleCellExperiment dyngen stringr dplyr
create_5_lineages_simu <- function(multiplier = 1, nSim = 200) {
  # nSim <- 20; multiplier <- 1
  set.seed(15)
  bblego_list <- list(
    bblego_start("A", type = "simple", num_modules = 2),
    bblego_linear("A", to = "B", num_modules = 2),
    bblego_branching("B", c("C", "D", "E")),
    bblego_branching("C", c("F", "G")),
    bblego_branching("E", c("H", "I")),
    bblego_end("D", num_modules = 4), 
    bblego_end("F", num_modules = 4), bblego_end("G", num_modules = 4),
    bblego_end("H", num_modules = 6), bblego_end("I", num_modules = 6)
  )
  backbone <- bblego(.list = bblego_list)
  model_common <- .generate_model_common(backbone, nSim)
  model_comb <- .generate_model_comb(model_common, "B17", multiplier, nSim) %>%
    .calculate_dimred()
  sce <- .convert_to_sce(model_comb)
  return(sce)
}

