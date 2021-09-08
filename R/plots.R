#' @title Plot simulation dataset by condition and cluster
#' @description Takes as input the low dimension representation and the
#' skeleton structure of WT and KO and plot them, colored by cluster
#' @param model Output of \code{\link{initialise_model}}
#' @export
#' @details Code mostly copied from https://github.com/dynverse/dyngen/blob/master/R/plotting.R
#' @import ggplot2 tidygraph ggraph dplyr grid stringr dynutils tibble
#' @importFrom igraph layout.graphopt
plot_backbone_modulenet_simplify <- function (model) {
  module_id <- color <- from <- to <- strength <- effect <- name <- NULL
  node_legend <- model$backbone$module_info %>%
    dplyr::select(module_id, color) %>% tibble::deframe()
  nodes <- model$backbone$module_info %>% dplyr::rename(name = module_id)
  edges <- model$backbone$module_network %>% arrange(from == to)
  edges <- edges %>% filter(to %in% c(edges$from, "D4", "C5")) %>%
    filter(!str_detect(from, "Burn"))
  nodes <- nodes %>% filter(name %in% c(edges$from, edges$to))
  gr <- tidygraph::tbl_graph(nodes = nodes, edges = edges)
  layout <- gr %>% igraph::layout.graphopt(charge = 0.01, niter = 10000) %>%
    dynutils::scale_minmax() %>% as.data.frame()
  rownames(layout) <- nodes$name
  colnames(layout) <- c("x", "y")
  r <- 0.03
  cap <- ggraph::circle(4, "mm")
  str <- 0.2
  arrow_up <- grid::arrow(type = "closed", angle = 30,
                          length = grid::unit(3, "mm"))
  arrow_down <- grid::arrow(type = "closed", angle = 89,
                            length = grid::unit(3, "mm"))
  p <- ggraph::ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
    ggraph::geom_edge_loop(aes(width = strength, strength = str,
                       filter = effect >= 0), arrow = arrow_up, start_cap = cap,
                   end_cap = cap) +
    ggraph::geom_edge_loop(aes(width = strength,
                       strength = str, filter = effect < 0), arrow = arrow_down,
                   start_cap = cap, end_cap = cap) +
    ggraph::geom_edge_fan(aes(width = strength,
                      filter = effect >= 0), arrow = arrow_up, start_cap = cap,
                  end_cap = cap) +
    ggraph::geom_edge_fan(aes(width = strength, filter = effect < 0),
                  arrow = arrow_down, start_cap = cap,
                  end_cap = cap) +
    ggraph::geom_node_circle(aes(r = r, colour = name),  fill = "white") +
    ggraph::geom_node_text(aes(label = name)) +
    ggraph::theme_graph(base_family = "Helvetica") +
    scale_colour_manual(values = node_legend) +
    ggraph::scale_edge_width_continuous(trans = "log10",
                                range = c(0.5, 3)) +
    coord_equal()
  return(p)
}

#' @title Plot simulation dataset by condition and cluster
#' @description Takes as input the low dimension representation and the
#' skeleton structure of WT and KO and plot them, colored by cluster
#' @param sce A single cell experiment object, result of the simulation
#' @export
#' @import ggplot2 SingleCellExperiment
#' @importFrom dplyr filter bind_cols
plot_reduced_dim <- function(sce) {
  df <- bind_cols(
    colData(sce) %>% as.data.frame(),
    reducedDim(sce) %>% as.data.frame()
  )
  p <- ggplot(df) +
    geom_point(size = .5, aes(x = comp_1, y = comp_2, col = from)) +
    guides(col = guide_legend(override.aes = list(size = 2))) +
    labs(col = "Cluster", x= "Dim1", y = "Dim2") +
    scale_x_continuous(limits = c(-.5, .5)) +
    theme_classic()  +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    facet_wrap(~condition)
  return(p)
}

#' @title Plot simulation dataset by condition
#' @description Takes as input the low dimension representation and the
#' skeleton structure of WT and KO and plot them, colored by condition.
#' @param sce A single cell experiment object, result of the simulation
#' @param sample Fraction of the data to plot. Default to 1.
#' @export
#' @import ggplot2 SingleCellExperiment
#' @importFrom dplyr filter bind_cols sample_frac
plot_reduced_dim_together <- function(sce, sample = 1) {
  colors <- c("#0072B2", "#D55E00", "#F0E442")
  df <- dplyr::bind_cols(colData(sce) %>% as.data.frame(),
                         reducedDim(sce) %>% as.data.frame()) %>%
    dplyr::group_by(condition) %>% 
    dplyr::sample_frac(size = sample) %>%
    dplyr::ungroup() %>%
    dplyr::sample_frac(size = 1)
  p <- ggplot(df, aes(col = condition)) +
    geom_point(size = .1, aes(x = comp_1, y = comp_2)) +
    scale_color_manual(values = colors, breaks = c("WT", "KO", "UP")) +
    labs(col = "Conditions", x = "Dim1", y = "Dim2") +
    scale_x_continuous(limits = c(-.5, .5)) +
    theme_classic()  +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  return(p)
}

#' @title Plot examples for figure 1
#' @description Takes as input the low dimension representation and the
#' skeleton structure of WT and KO and plot them.
#' @param sd_ko Low dim rep for KO
#' @param mst_ko network structure for KO
#' @param sd_wt Low dim rep for WT
#' @param mst_wt network structure for WT
#' @param title Plot title
#' @param colors color scheme
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
plot_example <- function(sd_ko, mst_ko, sd_wt, mst_wt, title, colors) {
  p <- ggplot(mapping = aes(x = Dim1, y = Dim2, col = conditions)) +
    geom_point(data = mst_wt, size = 4) +
    geom_path(data = mst_wt, mapping = aes(group = lineages), size = 1.5) +
    geom_point(data = mst_ko, size = 4) +
    geom_path(data = mst_ko, mapping = aes(group = lineages), size = 1.5) +
    geom_point(data = sd_wt, alpha = .4) +
    geom_point(data = sd_ko %>% dplyr::filter(conditions == "B"), alpha = .8) +
    theme_classic() +
    labs(col = "Condition") +
    scale_color_manual(values = colors) +
    guides(col = FALSE) +
    theme(legend.position = c(.8, .5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = .5, size = 15)) +
    ggtitle(title)
  return(p)
}
