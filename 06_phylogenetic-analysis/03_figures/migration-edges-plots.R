# migration edges plots
## plot optm migration edges results
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(OptM)) install.packages("OptM")
if (!require(patchwork)) install.packages("patchwork")

### function to read and plot data
optm_plots <- function(data = "hd"){
  # read treemix results and run optm analysis
  results <- OptM::optM(paste0("06_phylogenetic-analysis/02_output/01_migration-edges/",
                               dplyr::case_when(data == "hd" ~ "01_hd",
                                                data == "ld" ~ "02_ld"))) %>%
    dplyr::rename_with(~ gsub("\\(|)",
                              '',
                              .x))
  # plot Lm results
  a <- results %>%
    dplyr::mutate(lowerLm = meanLm - sdLm,
                  upperLm = meanLm + sdLm) %>%
    ggplot2::ggplot(ggplot2::aes(x = m,
                                 y = meanLm,
                                 ymin = lowerLm,
                                 ymax = upperLm)) +
    ggplot2::geom_pointrange() +
    ggplot2::labs(y = "Mean *L*(*m*)+/-SD") +
    ggplot2::scale_x_continuous(breaks = c(1:15),
                                expand = c(0,
                                           0),
                                limits = c(0.9,
                                           15.1)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggtext::element_markdown())
  # plot f results
  b <- results %>%
    dplyr::mutate(lowerf = meanf - sdf,
                  upperf = meanf + sdf) %>%
    ggplot2::ggplot(ggplot2::aes(x = m,
                                 y = meanf,
                                 ymin = lowerf,
                                 ymax = upperf)) +
    ggplot2::geom_hline(yintercept = 0.998,
                        linetype = 2) +
    ggplot2::geom_pointrange() +
    ggplot2::labs(y = "Variance explained") +
    ggplot2::scale_x_continuous(breaks = c(1:15),
                                expand = c(0,
                                           0),
                                limits = c(0.9,
                                           15.1)) +
    ggplot2::annotate("text",
                      x = 1.5,
                      y = 0.994,
                      label = "0.998",
                      size = 3.5) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  # plot deltam results
  c <- results %>%
    ggplot2::ggplot(ggplot2::aes(x = m,
                                 y = Deltam)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = "*m* (migration edges)",
                  y = expression(Delta*italic("m"))) +
    ggplot2::scale_x_continuous(breaks = c(1:15),
                                expand = c(0,
                                           0),
                                limits = c(0.9,
                                           15.1)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown())
  # combine plots
  plot <- patchwork::wrap_plots(list(a,
                                     b,
                                     c),
                                ncol = 1) +
    patchwork::plot_annotation(tag_levels = 'A')
  # save plots
  png(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-migration-edges.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-migration-edges.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### apply function to hd and ld data
optm_plots()
optm_plots(data = "ld")