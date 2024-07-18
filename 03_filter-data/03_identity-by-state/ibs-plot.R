# ibs plot
## plot ibs heatmap
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(khroma)) install.packages("khroma")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(viridis)) install.packages("viridis")

### read results
ibs <- readr::read_table("03_filter-data/03_identity-by-state/ars-mind.mibs",
                         col_names = F) %>%
  as.data.frame()

### read ids
ids <- readr::read_table("03_filter-data/03_identity-by-state/ars-mind.mibs.id",
                         col_names = F) %>%
  dplyr::mutate(population = X1,
                id = stringr::str_replace_all(X2,
                                              "\\.",
                                              "_"))

### set ids as column names
colnames(ibs) <- ids$id

### set ids as row names
rownames(ibs) <- ids$id

### set populations
populations <- c("HOLS", "ANGU", "JERS",
                 "ROMA", "CHIA", "MARC", "MARE", "ALEN",
                 "MUTU", "LAGU", "NDAG",
                 "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                 "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                 "THAR", "GIR", "NELO")

### set group labels
group_labels <- c("European<br>*Bos taurus*",
                  "European<br>hybrid",
                  "African<br>*Bos taurus*",
                  "Trypanotolerant<br>African hybrid",
                  "Trypanosusceptible<br>African hybrid",
                  "*Bos<br>indicus*")

### set group levels
group_levels <- c("European Bos taurus",
                  "European hybrid",
                  "African Bos taurus",
                  "Trypanotolerant African hybrid",
                  "Trypanosusceptible African hybrid",
                  "Bos indicus")

### sort ids
ids_sort <- ids %>%
  dplyr::transmute(population = factor(population,
                                       levels = populations),
                   id = id) %>%
  dplyr::arrange(population,
                 id)

### sort results by ids
ibs_sort <- ibs[match(ids_sort$id, rownames(ibs)), ] %>%
  dplyr::select(ids_sort$id)

### identify triangular matrix
tri <- lower.tri(ibs_sort,
                 diag = T)

### extract triangular matrix values
ibs_tri <- data.frame(id1 = colnames(ibs_sort)[col(ibs_sort)[tri]],
                      id2 = rownames(ibs_sort)[row(ibs_sort)[tri]],
                      value = ibs_sort[tri])

### rename id columns
ids_sort <- ids_sort %>%
  dplyr::rename(population1 = population,
                id1 = id)

### add population names for id1
ibs_tri <- ibs_tri %>%
  dplyr::inner_join(ids_sort)

### rename id columns
ids_sort <- ids_sort %>%
  dplyr::rename(population2 = population1,
                id2 = id1)

### add population names for id2
ibs_tri <- ibs_tri %>%
  dplyr::inner_join(ids_sort)

### calculate mean values
ibs_tri_mean <- ibs_tri %>%
  dplyr::mutate(population1 = factor(population1,
                                     levels = populations),
                population2 = factor(population2,
                                     levels = populations),
                group1 = dplyr::case_when(population1 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                          population1 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                          population1 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                          population1 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                          population1 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                          population1 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus"),
                group2 = dplyr::case_when(population2 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                          population2 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                          population2 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                          population2 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                          population2 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                          population2 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus"),
                group1 = factor(group1,
                                labels = group_labels,
                                levels = group_levels),
                group2 = factor(group2,
                                labels = group_labels,
                                levels = group_levels)) %>%
  dplyr::group_by(population1,
                  population2) %>%
  dplyr::mutate(mean = mean(value))

### set custom tick lengths
ticks <- dplyr::tibble(population1 = populations,
                       length = c(rep(c(0.3,
                                        -0.3),
                                      12))) %>%
  dplyr::mutate(group1 = dplyr::case_when(population1 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                          population1 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                          population1 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                          population1 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                          population1 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                          population1 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus"),
                population2 = population1,
                group2 = group1)

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### generate plot
plot <- ibs_tri_mean %>%
  ggplot2::ggplot(ggplot2::aes(x = interaction(population1,
                                               group1),
                               y = interaction(population2,
                                               group2),
                               fill = mean)) +
  ggplot2::geom_raster() +
  ggplot2::scale_x_discrete(expand = c(0,
                                       0),
                            guide = ggh4x::guide_axis_nested(n.dodge = 2)) +
  ggplot2::scale_y_discrete(expand = c(0,
                                       0),
                            guide = ggh4x::guide_axis_nested(),
                            limits = rev) +
  viridis::scale_fill_viridis(direction = -1,
                              limits = c(0,
                                         1),
                              name = "Mean\nIBS") +
  ggplot2::guides(fill = ggplot2::guide_colourbar(draw.ulim = F,
                                                  draw.llim = F)) +
  ggplot2::geom_linerange(data = ticks,
                          ggplot2::aes(x = interaction(factor(population1),
                                                       factor(group1,
                                                              labels = group_labels,
                                                              levels = group_levels)),
                                       y = interaction(factor(population2),
                                                       factor(group2,
                                                              labels = group_labels,
                                                              levels = group_levels)),
                                       ymax = length,
                                       ymin = 0.5),
                          colour = "grey70",
                          inherit.aes = F,
                          linewidth = ggplot2::rel(0.25)) +
  ggplot2::coord_equal(clip = "off",
                       expand = F,
                       ylim = c(0.5,
                                24.5)) +
  ggplot2::theme_light() +
  ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(),
                 ggh4x.axis.nestline.x = ggplot2::element_line(colour = c(pal[17],
                                                                          pal[15],
                                                                          pal[9],
                                                                          pal[5],
                                                                          pal[26],
                                                                          pal[28]),
                                                               linewidth = 0.5),
                 ggh4x.axis.nestline.y = ggplot2::element_line(colour = rev(c(pal[17],
                                                                              pal[15],
                                                                              pal[9],
                                                                              pal[5],
                                                                              pal[26],
                                                                              pal[28])),
                                                               linewidth = 0.5),
                 ggh4x.axis.nesttext.x = ggtext::element_markdown(colour = "black",
                                                                  size = ggplot2::rel(1.1)),
                 ggh4x.axis.nesttext.y = ggtext::element_markdown(colour = "black",
                                                                  size = ggplot2::rel(1.1),
                                                                  hjust = 0.5),
                 legend.background = ggplot2::element_blank(),
                 legend.position = c(0.935,
                                     0.84),
                 legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank())

### save plot
png("03_filter-data/03_identity-by-state/ibs.png",
    height = 6,
    res = 300,
    units = "in",
    width = 7)
print(plot)
dev.off()
pdf("03_filter-data/03_identity-by-state/ibs.pdf",
    height = 6,
    width = 7)
print(plot)
dev.off()