
.create_sticker_image <- function(loc) {
  c1 <- condiments::create_differential_topology(shift = 0, noise = 0)$sd
  c2 <- c1 %>%
    dplyr::mutate(Dim2 = Dim2 - .5)
  cs <- dplyr::bind_rows(c1 = c1, c2 = c2, .id = "cond")
  c1 <- condiments::create_differential_topology(shift = 0, n_cells = 2000)$sd
  c2 <- c1 %>%
    dplyr::mutate(Dim2 = Dim2 - .5)
  cs_pts <- dplyr::bind_rows(c1 = c1, c2 = c2, .id = "cond") %>%
    dplyr::sample_frac(1)

  p <- ggplot2::ggplot(cs, ggplot2::aes(x = Dim1, y = Dim2, col = cond)) +
    ggplot2::geom_point(size = 1, data = cs_pts, alpha = .7) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(lineages, cond)),
                       size = 5) +
    ggplot2::theme_void() +
    ggplot2::scale_color_manual(values = c("c1" = "#CE2522", "c2" = "#FEDC56")) +
    ggplot2::guides(col = FALSE)
  ggplot2::ggsave(filename = loc, plot = p,
                  bg = "transparent", width = 8, height = 6)
}
