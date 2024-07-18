# structure plots
## plot structure results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(khroma)) install.packages("khroma")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(tidyr)) install.packages("tidyr")

### function to read and plot structure data
structure_plot <- function(i, colours, data = "hd", order, tick1, tick2){
  # read and edit fam file
  fam <- readr::read_delim(paste0("03_filter-data/05_filter-snps/",
                                  data,
                                  ".fam"),
                           col_names = F,
                           show_col_types = F) %>%
    dplyr::rename(population = X1,
                  id = X2) %>%
    dplyr::mutate(group = dplyr::case_when(population %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                           population %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                           population %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                           population %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                           population %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                           population %in% c("THAR", "GIR", "NELO") ~ "Bos indicus")) %>%
    dplyr::mutate(population = factor(population,
                                      levels = c("HOLS", "ANGU", "JERS",
                                                 "ROMA", "CHIA", "MARC", "MARE", "ALEN",
                                                 "MUTU", "LAGU", "NDAG",
                                                 "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                                                 "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                                                 "THAR", "GIR", "NELO")),
                  group = factor(group,
                                 labels = c("European<br>*Bos taurus*",
                                            "European<br>hybrid",
                                            "African<br>*Bos taurus*",
                                            "Trypanotolerant<br>African hybrid",
                                            "Trypanosusceptible<br>African hybrid",
                                            "*Bos<br>indicus*"),
                                 levels = c("European Bos taurus",
                                            "European hybrid",
                                            "African Bos taurus",
                                            "Trypanotolerant African hybrid",
                                            "Trypanosusceptible African hybrid",
                                            "Bos indicus")))
  # set length of custom axis ticks
  ticks <- dplyr::tibble(population = fam$population,
                         id = fam$id,
                         group = fam$group,
                         length = c(rep(NA, 2), tick2, rep(NA, 3),
                                    rep(NA, 12), tick2, rep(NA, 13),
                                    rep(NA, 12), tick1, rep(NA, 12),
                                    rep(NA, 44), tick1, rep(NA, 45),
                                    rep(NA, 24), tick1, rep(NA, 25),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 15), tick1, rep(NA, 95),
                                    rep(NA, 13), tick1, rep(NA, 14),
                                    rep(NA, 17), tick1, rep(NA, 18),
                                    rep(NA, 11), tick1, rep(NA, 11),
                                    rep(NA, 7), tick2, rep(NA, 8),
                                    rep(NA, 10), tick1, rep(NA, 11),
                                    rep(NA, 2), tick2, rep(NA, 2),
                                    rep(NA, 6), tick2, rep(NA, 6),
                                    rep(NA, 2), tick1, rep(NA, 2),
                                    rep(NA, 9), tick1, rep(NA, 9),
                                    rep(NA, 23), tick1, rep(NA, 23),
                                    rep(NA, 20), tick2, rep(NA, 20),
                                    rep(NA, 18), tick2, rep(NA, 19),
                                    rep(NA, 13), tick2, rep(NA, 13),
                                    rep(NA, 25), tick2, rep(NA, 25),
                                    rep(NA, 7), tick2, rep(NA, 8),
                                    rep(NA, 11), tick2, rep(NA, 11),
                                    rep(NA, 6), tick2, rep(NA, 6)))
  # read and plot structure results
  readr::read_delim(paste0("05_structure-analysis/02_output/",
                           dplyr::case_when(data == "hd" ~ "01_hd",
                                            data == "ld" ~ "02_ld"),
                           "/fS_run_K.",
                           i,
                           ".meanQ"),
                    col_names = F,
                    delim = "  ",
                    show_col_types = F) %>%
    dplyr::mutate(population = fam$population,
                  id = fam$id,
                  group = fam$group) %>%
    dplyr::rename(all_of(order)) %>%
    tidyr::pivot_longer(cols = 1:all_of(i)) %>%
    ggplot2::ggplot(ggplot2::aes(x = id,
                                 y = value,
                                 fill = name)) +
    ggplot2::geom_bar(stat = "identity",
                      position = "fill",
                      width = 1) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0)) +
    ggh4x::facet_nested(~ group + population,
                        scales = "free_x",
                        space = "free_x",
                        switch = "x",
                        nest_line = ggplot2::element_line(colour = "grey70",
                                                          linewidth = ggplot2::rel(0.25)),
                        strip = ggh4x::strip_nested(clip = "off",
                                                    text_x = list(ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggtext::element_markdown(colour = "black",
                                                                                           size = ggplot2::rel(1.2)),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.6,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 1,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.6,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.3,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.8,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.4,
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        hjust = 0.3,
                                                                                        vjust = 0.2),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.8),
                                                                  ggplot2::element_text(colour = "grey30",
                                                                                        vjust = 0.2)))) +
    ggplot2::geom_linerange(data = ticks,
                            ggplot2::aes(x = id,
                                         ymax = length,
                                         ymin = 0),
                            colour = "grey70",
                            inherit.aes = F,
                            linewidth = ggplot2::rel(0.25),
                            na.rm = T) +
    ggplot2::coord_cartesian(clip = "off",
                             expand = F,
                             ylim = c(0,
                                      1)) +
    ggplot2::ggtitle(paste0("*K* = ",
                            i)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)),
                   plot.title = ggtext::element_markdown(size = ggplot2::rel(1)),
                   plot.margin = ggplot2::unit(c(5.5,
                                                 5.5,
                                                 -2.5,
                                                 5.5),
                                               "pt"),
                   strip.background = ggplot2::element_blank())
}

### function to save structure plot for individual k value
individual_k_plot <- function(i, colours, data = "hd", order, tick1 = -0.01, tick2 = -0.04){
  plot <- structure_plot(i,
                         colours = colours,
                         data = data,
                         order = order,
                         tick1 = tick1,
                         tick2 = tick2)
  png(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             data,
             "-k-",
             stringr::str_pad(i,
                              2,
                              pad = "0",
                              side = "left"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to draw structure plot for individual k value to combine with others
combined_k_plot <- function(i, colours, data = "hd", order, margin = 5.5, tick1 = 0, tick2 = 0){
  structure_plot(i,
                 colours = colours,
                 data = data,
                 order = order,
                 tick1 = tick1,
                 tick2 = tick2) +
    ggh4x::facet_nested(~ group + population,
                        scales = "free_x",
                        space = "free_x",
                        strip = ggh4x::strip_nested(text_x = list(ggplot2::element_blank()))) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0,
                                                 0,
                                                 margin,
                                                 0),
                                               "pt"))
}

### function to save combined structure plots
save_combined_plot <- function(k1 = 2, k2, data = "hd", height = 6){
  patch <- list()
  for (i in k1:(k2-1)){
    name <- paste0(data,
                   "_",
                   stringr::str_pad(i,
                                    2,
                                    side = "left",
                                    pad = "0"))
    patch[[name]] <- get(name)
  }
  len <- length(patch)
  patch[[len + 1]] <- get(paste0(data,
                                 "_",
                                 stringr::str_pad(k2,
                                                  2,
                                                  side = "left",
                                                  pad = "0"),
                                 "_axis"))
  plot <- patchwork::wrap_plots(patch,
                                ncol = 1) +
    patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::unit(c(5.5,
                                                                                    5.5,
                                                                                    -2.5,
                                                                                    5.5),
                                                                                  "pt")))
  png(paste0("05_structure-analysis/03_figures/",
             data,
             "-combined-k-",
             stringr::str_pad(k1,
                              2,
                              pad = "0",
                              side = "left"),
             "-",
             stringr::str_pad(k2,
                              2,
                              pad = "0",
                              side = "left"),
             ".png"),
      height = height,
      res = 600,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_structure-analysis/03_figures/",
             data,
             "-combined-k-",
             stringr::str_pad(k1,
                              2,
                              pad = "0",
                              side = "left"),
             "-",
             stringr::str_pad(k2,
                              2,
                              pad = "0",
                              side = "left"),
             ".pdf"),
      height = height,
      width = 9)
  print(plot)
  dev.off()
}

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### assign colours and orders
hd_02_colours <- c(pal[17], pal[28])
hd_02_order <- c(a = "X2", b = "X1")

ld_02_colours <- c(pal[17], pal[28])
ld_02_order <- c(a = "X2", b = "X1")

hd_03_colours <- c(pal[9], pal[17], pal[28])
hd_03_order <- c(a = "X1", b = "X3", c = "X2")

ld_03_colours <- c(pal[9], pal[17], pal[28])
ld_03_order <- c(a = "X1", b = "X2", c = "X3")

hd_04_colours <- c(pal[26], pal[9], pal[17], pal[28])
hd_04_order <- c(a = "X3", b = "X2", c = "X4", d = "X1")

ld_04_colours <- c(pal[26], pal[9], pal[17], pal[28])
ld_04_order <- c(a = "X3", b = "X4", c = "X2", d = "X1")

hd_05_colours <- c(pal[15], pal[26], pal[9], pal[17], pal[28])
hd_05_order <- c(a = "X3", b = "X2", c = "X1", d = "X5", e = "X4")

ld_05_colours <- c(pal[15], pal[26], pal[9], pal[17], pal[28])
ld_05_order <- c(a = "X1", b = "X2", c = "X3", d = "X5", e = "X4")

hd_06_colours <- c(pal[7], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
hd_06_order <- c(a = "X5", b = "X3", c = "X6", d = "X4", e = "X2",
                 f = "X1")

ld_06_colours <- c(pal[19], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
ld_06_order <- c(a = "X6", b = "X3", c = "X1", d = "X5", e = "X4",
                 f = "X2")

hd_07_colours <- c(pal[25], pal[19], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
hd_07_order <- c(a = "X6", b = "X2", c = "X3", d = "X7", e = "X1",
                 f = "X5", g = "X4")

ld_07_colours <- c(pal[25], pal[19], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
ld_07_order <- c(a = "X2", b = "X4", c = "X1", d = "X6", e = "X3",
                 f = "X7", g = "X5")

hd_08_colours <- c(pal[5], pal[25], pal[7], pal[15], pal[26],
                   pal[9], pal[17], pal[28])
hd_08_order <- c(a = "X1", b = "X7", c = "X8", d = "X6", e = "X4",
                 f = "X2", g = "X5", h = "X3")

ld_08_colours <- c(pal[7], pal[25], pal[19], pal[15], pal[26],
                   pal[9], pal[17], pal[28])
ld_08_order <- c(a = "X5", b = "X7", c = "X2", d = "X6", e = "X3",
                 f = "X4", g = "X1", h = "X8")

hd_09_colours <- c(pal[19], pal[5], pal[25], pal[7], pal[15],
                   pal[26], pal[9], pal[17], pal[28])
hd_09_order <- c(a = "X3", b = "X5", c = "X2", d = "X1", e = "X9",
                 f = "X8", g = "X6", h = "X7", i = "X4")

ld_09_colours <- c(pal[18], pal[30], pal[25], pal[7], pal[15],
                   pal[26], pal[9], pal[17], pal[28])
ld_09_order <- c(a = "X2", b = "X9", c = "X1", d = "X7", e = "X4",
                 f = "X5", g = "X3", h = "X6", i = "X8")

hd_10_colours <- c(pal[14], pal[19], pal[4], pal[25], pal[7],
                   pal[15], pal[26], pal[9], pal[17], pal[28])
hd_10_order <- c(a = "X6", b = "X8", c = "X2", d = "X1", e = "X4",
                 f = "X3", g = "X9", h = "X7", i = "X10", j = "X5")

ld_10_colours <- c(pal[18], pal[19], pal[5], pal[25], pal[7],
                   pal[15], pal[26], pal[9], pal[17], pal[28])
ld_10_order <- c(a = "X4", b = "X5", c = "X1", d = "X8", e = "X10",
                 f = "X9", g = "X2", h = "X3", i = "X7", j = "X6")

hd_11_colours <- c(pal[30], pal[14], pal[19], pal[4], pal[25],
                   pal[7], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
hd_11_order <- c(a = "X9", b = "X4", c = "X5", d = "X2", e = "X8",
                 f = "X1", g = "X7", h = "X10", i = "X6", j = "X3",
                 k = "X11")

ld_11_colours <- c(pal[4], pal[14], pal[19], pal[5], pal[25],
                   pal[7], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
ld_11_order <- c(a = "X8", b = "X1", c = "X7", d = "X3", e = "X10",
                 f = "X4", g = "X11", h = "X6", i = "X2", j = "X9",
                 k = "X5")

hd_12_colours <- c(pal[24], pal[4], pal[14], pal[19], pal[5],
                   pal[25], pal[7], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
hd_12_order <- c(a = "X5", b = "X3", c = "X12", d = "X8", e = "X9",
                 f = "X7", g = "X11", h = "X6", i = "X4", j = "X10",
                 k = "X1", l = "X2")

ld_12_colours <- c(pal[23], pal[30], pal[14], pal[19], pal[5],
                   pal[25], pal[7], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
ld_12_order <- c(a = "X12", b = "X1", c = "X5", d = "X8", e = "X4",
                 f = "X7", g = "X11", h = "X10", i = "X6", j = "X9",
                 k = "X3", l = "X2")

hd_13_colours <- c(pal[1], pal[18], pal[30], pal[14], pal[19],
                   pal[5], pal[25], pal[7], pal[15], pal[26],
                   pal[9], pal[17], pal[28])
hd_13_order <- c(a = "X13", b = "X10", c = "X5", d = "X7", e = "X3",
                 f = "X9", g = "X8", h = "X6", i = "X4", j = "X2",
                 k = "X1", l = "X11", m = "X12")

ld_13_colours <- c(pal[22], pal[24], pal[18], pal[4], pal[14],
                   pal[19], pal[25], pal[7], pal[15], pal[26],
                   pal[9], pal[17], pal[28])
ld_13_order <- c(a = "X6", b = "X13", c = "X4", d = "X11", e = "X5",
                 f = "X1", g = "X10", h = "X8", i = "X3", j = "X12",
                 k = "X9", l = "X7", m = "X2")

hd_14_colours <- c(pal[29], pal[23], pal[30], pal[4], pal[14],
                   pal[19], pal[5], pal[25], pal[7], pal[15],
                   pal[26], pal[9], pal[17], pal[28])
hd_14_order <- c(a = "X4", b = "X10", c = "X2", d = "X9", e = "X1",
                 f = "X5", g = "X12", h = "X7", i = "X8", j = "X6",
                 k = "X14", l = "X11", m = "X13", n = "X3")

ld_14_colours <- c(pal[24], pal[23], pal[18], pal[4], pal[30],
                   pal[19], pal[5], pal[25], pal[7], pal[15],
                   pal[26], pal[9], pal[17], pal[28])
ld_14_order <- c(a = "X10", b = "X12", c = "X9", d = "X7", e = "X5",
                 f = "X13", g = "X6", h = "X2", i = "X11", j = "X14",
                 k = "X1", l = "X8", m = "X3", n = "X4")

hd_15_colours <- c(pal[23], pal[24], pal[1], pal[18], pal[30],
                   pal[4], pal[14], pal[19], pal[25], pal[7],
                   pal[15], pal[26], pal[9], pal[17], pal[28])
hd_15_order <- c(a = "X3", b = "X13", c = "X15", d = "X10", e = "X2",
                 f = "X12", g = "X4", h = "X14", i = "X11", j = "X6",
                 k = "X5", l = "X8", m = "X9", n = "X7", o = "X1")

ld_15_colours <- c(pal[24], pal[22],  pal[18], pal[4], pal[30],
                   pal[14], pal[19], pal[5], pal[25], pal[7],
                   pal[15], pal[26], pal[9], pal[17], pal[28])
ld_15_order <- c(a = "X13", b = "X11", c = "X3", d = "X8", e = "X9",
                 f = "X7", g = "X5", h = "X14", i = "X6", j = "X4",
                 k = "X10", l = "X1", m = "X12", n = "X15", o = "X2")

hd_16_colours <- c(pal[13], pal[23], pal[24], pal[18], pal[30],
                   pal[4], pal[14], pal[19], pal[5], pal[25],
                   pal[7], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
hd_16_order <- c(a = "X9", b = "X7", c = "X4", d = "X2", e = "X11",
                 f = "X1", g = "X8", h = "X5", i = "X13", j = "X14",
                 k = "X6", l = "X12", m = "X10", n = "X16", o = "X3",
                 p = "X15")

ld_16_colours <- c(pal[11], pal[23], pal[24], pal[18], pal[4],
                   pal[30], pal[14], pal[19], pal[5], pal[25],
                   pal[7], pal[15], pal[26], pal[9], pal[17],
                   pal[28])
ld_16_order <- c(a = "X12", b = "X1", c = "X9", d = "X6", e = "X2",
                 f = "X5", g = "X7", h = "X11", i = "X15", j = "X8",
                 k = "X16", l = "X14", m = "X10", n = "X4", o = "X3",
                 p = "X13")

hd_17_colours <- c(pal[22], pal[13], pal[24], pal[1], pal[18],
                   pal[30], pal[4], pal[14], pal[19], pal[5],
                   pal[25], pal[7], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
hd_17_order <- c(a = "X2", b = "X14", c = "X1", d = "X13", e = "X4",
                 f = "X10", g = "X12", h = "X17", i = "X3", j = "X15",
                 k = "X16", l = "X9", m = "X7", n = "X5", o = "X6",
                 p = "X8", q = "X11")

ld_17_colours <- c("black", pal[29], pal[23], pal[24], pal[18],
                   pal[4], pal[30], pal[14], pal[19], pal[5],
                   pal[25], pal[7], pal[15], pal[26], pal[9],
                   pal[17], pal[28])
ld_17_order <- c(a = "X16", b = "X12", c = "X2", d = "X10", e = "X5",
                 f = "X9", g = "X7", h = "X4", i = "X1", j = "X14",
                 k = "X15", l = "X17", m = "X13", n = "X3", o = "X6",
                 p = "X8", q = "X11")

### apply functions to results
#### save individual results
##### k = 2
individual_k_plot(2,
                  colours = hd_02_colours,
                  order = hd_02_order)

individual_k_plot(2,
                  colours = ld_02_colours,
                  data = "ld",
                  order = ld_02_order)

##### k = 3
individual_k_plot(3,
                  colours = hd_03_colours,
                  order = hd_03_order)

individual_k_plot(3,
                  colours = ld_03_colours,
                  data = "ld",
                  order = ld_03_order)

##### k = 4
individual_k_plot(4,
                  colours = hd_04_colours,
                  order = hd_04_order)

individual_k_plot(4,
                  colours = ld_04_colours,
                  data = "ld",
                  order = ld_04_order)

##### k = 5
individual_k_plot(5,
                  colours = hd_05_colours,
                  order = hd_05_order)

individual_k_plot(5,
                  colours = ld_05_colours,
                  data = "ld",
                  order = ld_05_order)

##### k = 6
individual_k_plot(6,
                  colours = hd_06_colours,
                  order = hd_06_order)

individual_k_plot(6,
                  colours = ld_06_colours,
                  data = "ld",
                  order = ld_06_order)

##### k = 7
individual_k_plot(7,
                  colours = hd_07_colours,
                  order = hd_07_order)

individual_k_plot(7,
                  colours = ld_07_colours,
                  data = "ld",
                  order = ld_07_order)

##### k = 8
individual_k_plot(8,
                  colours = hd_08_colours,
                  order = hd_08_order)

individual_k_plot(8,
                  colours = ld_08_colours,
                  data = "ld",
                  order = ld_08_order)

##### k = 9
individual_k_plot(9,
                  colours = hd_09_colours,
                  order = hd_09_order)

individual_k_plot(9,
                  colours = ld_09_colours,
                  data = "ld",
                  order = ld_09_order)

##### k = 10
individual_k_plot(10,
                  colours = hd_10_colours,
                  order = hd_10_order)

individual_k_plot(10,
                  colours = ld_10_colours,
                  data = "ld",
                  order = ld_10_order)

##### k = 11
individual_k_plot(11,
                  colours = hd_11_colours,
                  order = hd_11_order)

individual_k_plot(11,
                  colours = ld_11_colours,
                  data = "ld",
                  order = ld_11_order)

##### k = 12
individual_k_plot(12,
                  colours = hd_12_colours,
                  order = hd_12_order)

individual_k_plot(12,
                  colours = ld_12_colours,
                  data = "ld",
                  order = ld_12_order)

##### k = 13
individual_k_plot(13,
                  colours = hd_13_colours,
                  order = hd_13_order)

individual_k_plot(13,
                  colours = ld_13_colours,
                  data = "ld",
                  order = ld_13_order)

##### k = 14
individual_k_plot(14,
                  colours = hd_14_colours,
                  order = hd_14_order)

individual_k_plot(14,
                  colours = ld_14_colours,
                  data = "ld",
                  order = ld_14_order)

##### k = 15
individual_k_plot(15,
                  colours = hd_15_colours,
                  order = hd_15_order)

individual_k_plot(15,
                  colours = ld_15_colours,
                  data = "ld",
                  order = ld_15_order)

##### k = 16
individual_k_plot(16,
                  colours = hd_16_colours,
                  order = hd_16_order)

individual_k_plot(16,
                  colours = ld_16_colours,
                  data = "ld",
                  order = ld_16_order)

##### k = 17
individual_k_plot(17,
                  colours = hd_17_colours,
                  order = hd_17_order)

individual_k_plot(17,
                  colours = ld_17_colours,
                  data = "ld",
                  order = ld_17_order)

#### save combined results
##### k = 2 - 3
hd_02 <- combined_k_plot(2,
                         colours = hd_02_colours,
                         margin = 0.5,
                         order = hd_02_order)

ld_02 <- combined_k_plot(2,
                         data = "ld",
                         colours = ld_02_colours,
                         margin = 0.5,
                         order = ld_02_order)

hd_03_axis <- structure_plot(3,
                             colours = hd_03_colours,
                             order = hd_03_order,
                             tick1 = -0.02,
                             tick2 = -0.09)

ld_03_axis <- structure_plot(3,
                             colours = ld_03_colours,
                             data = "ld",
                             order = ld_03_order,
                             tick1 = -0.02,
                             tick2 = -0.09)

save_combined_plot(k2 = 3)

save_combined_plot(k2 = 3, data = "ld")

##### k = 2 - 9
hd_02 <- combined_k_plot(2,
                         colours = hd_02_colours,
                         order = hd_02_order)

ld_02 <- combined_k_plot(2,
                         colours = ld_02_colours,
                         data = "ld",
                         order = ld_02_order)

hd_03 <- combined_k_plot(3,
                         colours = hd_03_colours,
                         order = hd_03_order)

ld_03 <- combined_k_plot(3,
                         colours = ld_03_colours,
                         data = "ld",
                         order = ld_03_order)

hd_04 <- combined_k_plot(4,
                         colours = hd_04_colours,
                         order = hd_04_order)

ld_04 <- combined_k_plot(4,
                         colours = ld_04_colours,
                         data = "ld",
                         order = ld_04_order)

hd_05 <- combined_k_plot(5,
                         colours = hd_05_colours,
                         order = hd_05_order)

ld_05 <- combined_k_plot(5,
                         colours = ld_05_colours,
                         data = "ld",
                         order = ld_05_order)

hd_06 <- combined_k_plot(6,
                         colours = hd_06_colours,
                         order = hd_06_order)

ld_06 <- combined_k_plot(6,
                         colours = ld_06_colours,
                         data = "ld",
                         order = ld_06_order)

hd_07 <- combined_k_plot(7,
                         colours = hd_07_colours,
                         order = hd_07_order)

ld_07 <- combined_k_plot(7,
                         colours = ld_07_colours,
                         data = "ld",
                         order = ld_07_order)

hd_08 <- combined_k_plot(8,
                         colours = hd_08_colours,
                         order = hd_08_order)

ld_08 <- combined_k_plot(8,
                         colours = ld_08_colours,
                         data = "ld",
                         order = ld_08_order)

hd_09_axis <- structure_plot(9,
                             colours = hd_09_colours,
                             order = hd_09_order,
                             tick1 = -0.03,
                             tick2 = -0.16) +
  ggplot2::theme(plot.margin = ggplot2::unit(c(0,
                                               0,
                                               0,
                                               0),
                                             "pt"))

ld_09_axis <- structure_plot(9,
                             colours = ld_09_colours,
                             data = "ld",
                             order = ld_09_order,
                             tick1 = -0.03,
                             tick2 = -0.16) +
  ggplot2::theme(plot.margin = ggplot2::unit(c(0,
                                               0,
                                               0,
                                               0),
                                             "pt"))

save_combined_plot(k2 = 9,
                   height = 13.5)

save_combined_plot(k2 = 9,
                   data = "ld",
                   height = 13.5)

##### k = 10 - 17
hd_10 <- combined_k_plot(10,
                         colours = hd_10_colours,
                         order = hd_10_order)

ld_10 <- combined_k_plot(10,
                         colours = ld_10_colours,
                         data = "ld",
                         order = ld_10_order)

hd_11 <- combined_k_plot(11,
                         colours = hd_11_colours,
                         order = hd_11_order)

ld_11 <- combined_k_plot(11,
                         colours = ld_11_colours,
                         data = "ld",
                         order = ld_11_order)

hd_12 <- combined_k_plot(12,
                         colours = hd_12_colours,
                         order = hd_12_order)

ld_12 <- combined_k_plot(12,
                         colours = ld_12_colours,
                         data = "ld",
                         order = ld_12_order)

hd_13 <- combined_k_plot(13,
                         colours = hd_13_colours,
                         order = hd_13_order)

ld_13 <- combined_k_plot(13,
                         colours = ld_13_colours,
                         data = "ld",
                         order = ld_13_order)

hd_14 <- combined_k_plot(14,
                         colours = hd_14_colours,
                         order = hd_14_order)

ld_14 <- combined_k_plot(14,
                         colours = ld_14_colours,
                         data = "ld",
                         order = ld_14_order)

hd_15 <- combined_k_plot(15,
                         colours = hd_15_colours,
                         order = hd_15_order)

ld_15 <- combined_k_plot(15,
                         colours = ld_15_colours,
                         data = "ld",
                         order = ld_15_order)

hd_16 <- combined_k_plot(16,
                         colours = hd_16_colours,
                         order = hd_16_order)

ld_16 <- combined_k_plot(16,
                         colours = ld_16_colours,
                         data = "ld",
                         order = ld_16_order)

hd_17_axis <- structure_plot(17,
                             colours = hd_17_colours,
                             order = hd_17_order,
                             tick1 = -0.03,
                             tick2 = -0.16) +
  ggplot2::theme(plot.margin = ggplot2::unit(c(0,
                                               0,
                                               0,
                                               0),
                                             "pt"))

ld_17_axis <- structure_plot(17,
                             colours = ld_17_colours,
                             data = "ld",
                             order = ld_17_order,
                             tick1 = -0.03,
                             tick2 = -0.16) +
  ggplot2::theme(plot.margin = ggplot2::unit(c(0,
                                               0,
                                               0,
                                               0),
                                             "pt"))

save_combined_plot(k1 = 10,
                   k2 = 17,
                   height = 13.5)

save_combined_plot(k1 = 10,
                   k2 = 17,
                   data = "ld",
                   height = 13.5)