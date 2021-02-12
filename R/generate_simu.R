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
