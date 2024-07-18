# inbreeding outliers plot
## plot inbreeding results with outliers
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(khroma)) install.packages("khroma")
if (!require(readr)) install.packages("readr")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### set colour palette
cols <-  c(pal[17], pal[18], pal[19],
           pal[15], pal[14], pal[13], pal[12], pal[11],
           pal[9], pal[8], pal[7],
           pal[5], pal[4], pal[3], pal[2], pal[1],
           pal[22], pal[23], pal[24], pal[25], pal[26],
           pal[28], pal[29], pal[30])

### set limits to order colours
lims <- c("HOLS", "ANGU", "JERS",
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

### set custom tick lengths
ticks <- dplyr::tibble(FID = lims,
                       length = rep(c(-0.66,
                                      -0.69),
                                    12)) %>%
  dplyr::mutate(group = dplyr::case_when(FID %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                         FID %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                         FID %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                         FID %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                         FID %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                         FID %in% c("THAR", "GIR", "NELO") ~ "Bos indicus"))

### read and plot inbreeding data
plot <- readr::read_table("03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs.het") %>%
  dplyr::group_by(FID) %>%
  dplyr::mutate(group = dplyr::case_when(FID %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                         FID %in% c("ALEN", "CHIA", "MARC", "MARE", "ROMA") ~ "European hybrid",
                                         FID %in% c("LAGU", "NDAG", "MUTU") ~ "African Bos taurus",
                                         FID %in% c("BORG", "KETE", "NDAM", "SHEK", "SOMB") ~ "Trypanotolerant African hybrid",
                                         FID %in% c("ANKO", "BORA", "EASZ", "KARA", "NGAN") ~ "Trypanosusceptible African hybrid",
                                         FID %in% c("GIR", "NELO", "THAR") ~ "Bos indicus")) %>%
  ggplot2::ggplot(ggplot2::aes(interaction(factor(FID,
                                                  levels = lims),
                                           factor(group,
                                                  labels = group_labels,
                                                  levels = group_levels)),
                               F)) +
  ggplot2::geom_boxplot(fill = NA,
                        outlier.size = 2.5) +
  ggplot2::geom_dotplot(ggplot2::aes(colour = FID,
                                     fill = FID),
                        binaxis = "y",
                        binwidth = 0.003,
                        dotsize = 5,
                        stackdir = "center",
                        stackratio = 0.3) +
  ggplot2::geom_boxplot(fill = NA,
                        outlier.shape = NA) +
  ggplot2::scale_colour_manual(limits = lims,
                               values = cols) +
  ggplot2::scale_fill_manual(limits = lims,
                             values = cols) +
  ggplot2::scale_x_discrete(guide = ggh4x::guide_axis_nested(n.dodge = 2)) +
  ggplot2::geom_linerange(data = ticks,
                          ggplot2::aes(interaction(factor(FID),
                                                   factor(group,
                                                          labels = group_labels,
                                                          levels = group_levels)),
                                       ymax = length,
                                       ymin = -0.65),
                          colour = "grey70",
                          inherit.aes = F,
                          linewidth = ggplot2::rel(0.25)) +
  ggplot2::coord_cartesian(clip = "off",
                           expand = F,
                           ylim = c(-0.65,
                                    0.65)) +
  ggplot2::scale_y_continuous(breaks = c(-0.6,
                                         -0.3,
                                         0,
                                         0.3,
                                         0.6)) +
  ggplot2::theme_light() +
  ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_text(face = "italic"),
                 ggh4x.axis.nestline = ggplot2::element_line(colour = c(pal[17],
                                                                        pal[15],
                                                                        pal[9],
                                                                        pal[5],
                                                                        pal[26],
                                                                        pal[28]),
                                                             linewidth = 0.5),
                 ggh4x.axis.nesttext.x = ggtext::element_markdown(colour = "black",
                                                                  size = ggplot2::rel(1.2)),
                 legend.position = "none",
                 panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)),
                 strip.background = ggplot2::element_blank(),
                 strip.placement = "outside")

### save plot
png("03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()