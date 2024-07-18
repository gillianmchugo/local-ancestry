# pca plots
## plot pca results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(khroma)) install.packages("khroma")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(tidyr)) install.packages("tidyr")

### function to read and plot data
pca_eigen_plot <- function(pc2 = 2, pc1 = 1, data = "hd"){
  # read and plot eigenvalues from eval file
  eigen <- readr::read_table(paste0("04_principal-component-analysis/02_output/",
                                    data,
                                    ".eval"),
                             col_names = F) %>%
    dplyr::slice(1:10) %>%
    dplyr::mutate(prop = X1/sum(X1),
                  pc = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(pc,
                                 prop,
                                 fill = factor(pc))) +
    ggplot2::geom_bar(position = "dodge",
                      stat="identity") +
    ggplot2::scale_x_continuous(breaks = c(1:10),
                                expand = c(0,
                                           0),
                                limits = c(0.55,
                                           10.45)) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                limits = c(0,
                                           0.55)) +
    ggplot2::scale_fill_manual(values = c(sapply(1:10,
                                                 function(i)
                                                   dplyr::case_when(pc1 == i ~ "#000000",
                                                                    pc2 == i ~ "#000000",
                                                                    .default = "#666666")))) +
    ggplot2::labs(x = "Principal Component",
                  y = "Proportion of Variance") +
    ggplot2::coord_equal(ratio = (10.45-0.55)/0.55) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "none")
  # read pca from evec file, add blank rows for group names in legend
  evec <- readr::read_table(paste0("04_principal-component-analysis/02_output/",
                                   data,
                                   ".evec"),
                            col_names = F) %>%
    dplyr::slice(-1) %>%
    tidyr::separate_wider_delim(X1, ":",
                                names = c("group",
                                          "id")) %>%
    dplyr::add_row(group = c("European Bos taurus",
                             "European hybrid",
                             "African Bos taurus",
                             "Trypanotolerant African hybrid",
                             "Trypanosusceptible African hybrid",
                             "Bos indicus"))
  # set limits to order legend
  lims <- c("European Bos taurus",
            "HOLS", "ANGU", "JERS",
            "European hybrid",
            "ROMA", "CHIA", "MARC", "MARE", "ALEN",
            "African Bos taurus",
            "MUTU", "LAGU", "NDAG",
            "Trypanotolerant African hybrid",
            "NDAM", "BORG", "SOMB", "KETE", "SHEK",
            "Trypanosusceptible African hybrid",
            "ANKO", "NGAN", "EASZ", "KARA", "BORA",
            "Bos indicus",
            "THAR", "GIR", "NELO")
  # format labels for legend
  labs <- c("<span style='color:black'>European<br>*Bos taurus*</span>",
            "HOLS", "ANGU", "JERS",
            "<span style='color:black'>European<br>hybrid</span>",
            "ROMA", "CHIA", "MARC", "MARE", "ALEN",
            "<span style='color:black'>African<br>*Bos taurus*</span>",
            "MUTU", "LAGU", "NDAG",
            "<span style='color:black'>Trypanotolerant<br>African hybrid</span>",
            "NDAM", "BORG", "SOMB", "KETE", "SHEK",
            "<span style='color:black'>Trypanosusceptible<br>African hybrid</span>",
            "ANKO", "NGAN", "EASZ", "KARA", "BORA",
            "<span style='color:black'>*Bos<br>indicus*</span>",
            "THAR", "GIR", "NELO")
  # set point shapes
  pchs <- c(NA,
            21:23,
            NA,
            21:25,
            NA,
            21:23,
            NA,
            21:25,
            NA,
            21:25,
            NA,
            21:23)
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # set colour palette
  cols <-  c(NA,
             pal[17], pal[18], pal[19],
             NA,
             pal[15], pal[14], pal[13], pal[12], pal[11],
             NA,
             pal[9], pal[8], pal[7],
             NA,
             pal[5], pal[4], pal[3], pal[2], pal[1],
             NA,
             pal[22], pal[23], pal[24], pal[25], pal[26],
             NA,
             pal[28], pal[29], pal[30])
  # plot pca
  pca <- evec %>%
    ggplot2::ggplot(ggplot2::aes(.[[as.numeric(pc1)+2]],
                                 .[[as.numeric(pc2)+2]],
                                 col = group,
                                 fill = group,
                                 pch = group)) +
    ggplot2::geom_point(size = 2.5,
                        stroke = 0.5) +
    ggplot2::scale_colour_manual(labels = labs,
                                 limits = lims,
                                 values = cols) +
    ggplot2::scale_fill_manual(labels = labs,
                               limits = lims,
                               values = ggplot2::alpha(cols,
                                                       0.5)) +
    ggplot2::scale_shape_manual(labels = labs,
                                limits = lims,
                                values = pchs) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3),
                    fill = ggplot2::guide_legend(ncol = 3),
                    shape = ggplot2::guide_legend(ncol = 3)) +
    ggplot2::labs(x = paste0("Principal Component ",
                             pc1,
                             " (",
                             round(eigen$data$prop[pc1]*100,
                                   2),
                             "%)"),
                  y = paste0("Principal Component ",
                             pc2,
                             " (",
                             round(eigen$data$prop[pc2]*100,
                                   2),
                             "%)")) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.key.height = ggplot2::unit(0.5,
                                                     "cm"),
                   legend.key.width = ggplot2::unit(0.1,
                                                    "cm"),
                   legend.text = ggtext::element_markdown(colour = "grey30",
                                                          size = ggplot2::rel(0.65)),
                   legend.title = ggplot2::element_blank())
  # set layout
  layout <- "AAB
             AAC"
  # combine plots
  plot <- pca +
    patchwork::guide_area() +
    eigen +
    patchwork::plot_annotation(tag_levels = 'A') +
    patchwork::plot_layout(design = layout,
                           guides = "collect")
  plot
  # save plot
  png(paste0("04_principal-component-analysis/03_figures/",
             data,
             "-pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("04_principal-component-analysis/03_figures/",
             data,
             "-pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### apply function to hd and ld data for pc 2 to 10
lapply(2:10,
       pca_eigen_plot)
lapply(2:10,
       pca_eigen_plot,
       data = "ld")