.create_sticker <- function(loc) {
  # ketchup <- system.file("figures/ketchup.png", package="condimentsPaper")
  ketchup <- "inst/figures/ketchup.png"
  ketchup <- magick::image_read(ketchup) %>%
    magick::image_background("transparent") %>%
    magick::image_rotate(60)
  # mustard <- system.file("figures/mustard.png", package="condimentsPaper")
  mustard <- "inst/figures/mustard.png"
  mustard <- magick::image_read(mustard) %>%
    magick::image_background("transparent") %>%
    magick::image_rotate(240)
  c1 <- condiments::create_differential_topology(shift = 0, noise = 0)$sd
  c2 <- c1 %>%
    dplyr::mutate(Dim2 = Dim2 - .5)
  cs <- dplyr::bind_rows(c1 = c1, c2 = c2, .id = "cond")
  c1 <- condiments::create_differential_topology(shift = 0, n_cells = 2000)$sd
  c2 <- c1 %>%
    dplyr::mutate(Dim2 = Dim2 - .5)
  cs_pts <- dplyr::bind_rows(c1 = c1, c2 = c2, .id = "cond") %>%
    dplyr::sample_frac(1)

  p <- ggplot(cs, aes(x = Dim1, y = Dim2, col = cond)) +
    geom_point(size = 1, data = cs_pts, alpha = .7) +
    geom_line(aes(group = interaction(lineages, cond)),
              size = 5) +
    theme_void() +
    scale_color_manual(values = c("c1" = "#CE2522", "c2" = "#FEDC56")) +
    guides(col = FALSE)
  p1 <- magick::image_graph()
  print(p)
  dev.off()
  p1 <- p1 %>% magick::image_rotate(60)
  p2 <- ggdraw() +
    draw_image(p1, scale = 1.15) +
    draw_image(ketchup, width = .25, height = .25, x = .7, y = .1) +
    draw_image(mustard, width = .25, height = .25, x = .45, y = -.08)
  p3 <- magick::image_graph()
  print(p2)
  dev.off()
  p3 <- p3 %>%
    magick::image_background("transparent")
  sticker(p3, package = "condiments",
          p_size = 20, s_x = 1, s_y = 1.1, s_width = 2.2, s_height = 2.2,
          h_fill = "white", p_color = "black", h_color = "black",
          filename = loc)
}
