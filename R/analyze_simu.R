.discovery_rates <- function(df, cutoff) {
  df %>%
    group_by(Type, test_type) %>%
    summarise(level = mean(adjusted_p_value <= cutoff), .groups = "drop") %>%
    mutate(Type = if_else(Type == "Null", "TNR", "TPR"),
           level = if_else(Type == "TNR", 1 - level, level)) %>%
    pivot_wider(names_from = Type, values_from = level) %>%
    return()
}
.predictive_rates <- function(df, cutoff) {
  df %>%
    mutate(called = if_else(adjusted_p_value <= cutoff, "PPV", "NPV")) %>%
    group_by(called, test_type) %>%
    summarise(level = mean(Type == "Real"), .groups = "drop") %>%
    mutate(level = if_else(called == "NPV", 1 - level, level)) %>%
    pivot_wider(names_from = called, values_from = level) %>%
    return()
}

#' @title Compute all 5 metrics for a given FDR threshold
#' @param df a data.frame, the output of the simulations
#' @param cutoff The FDR nominal level to control
#' @export
#' @import stringr dplyr tidyr
all_metrics <- function(df, cutoff) {
  df <- full_join(.discovery_rates(df, cutoff),
                  .predictive_rates(df, cutoff)) %>%
    mutate(F1 = 2 * (PPV * TPR) / (PPV + TPR))
  if (!"DAseq" %in% df$test_type) {
    df <- full_join(df, data.frame(test_type = "DAseq"))
  }
  df <- df %>%
    pivot_longer(-test_type, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>%
    mutate(ranks = base::rank(-value, ties.method = "min"),
           scaled = value - max(value, na.rm = TRUE)) %>%
    ungroup()
  return(df)
}

#' @title Clean the output from the various simulation
#' @param df a data.frame, the output of the simulations
#' @param nSim number of simulations for each condition
#' @export
#' @import stringr dplyr tidyr
clean_results <- function(df) {
  df %>%
    arrange(multiplier) %>%
    select(-1) %>%
    mutate(test_type = case_when(is.na(test_type) ~ "milo",
                                 str_detect(test_type, "region") ~ "DAseq",
                                 TRUE ~ paste0("condiments_", test_type)),
           adjusted_p_value = case_when(#!is.na(SpatialFDR) ~ SpatialFDR,
             !is.na(FDR) ~ FDR,
             TRUE ~ p.value)) %>%
    select(test_type, multiplier, adjusted_p_value) %>%
    return()
}
