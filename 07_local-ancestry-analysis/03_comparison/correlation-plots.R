# correlation plots
## correlation plots of local ancestry results 
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(parallel)) install.packages("parallel")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")
if (!require(viridis)) install.packages("viridis")

### function to draw correlation plots
correlation_plot <- function(analysis1, analysis2 = analysis1, software1 = "mosaic", software2 = "elai", data1 = "hd", data2 = data1){
  # read input data 1
  input1 <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                   dplyr::case_when(software1 == "mosaic" ~ "01_mosaic/",
                                                    software1 == "elai" ~ "02_elai/"),
                                   "02_output/",
                                   dplyr::case_when(data1 == "hd" ~ "01_hd/",
                                                    data1 == "ld" ~ "02_ld/"),
                                   dplyr::case_when(analysis1 == "ALEN" ~ "01_ALEN/",
                                                    analysis1 == "ANKO" ~ "02_ANKO/",
                                                    analysis1 == "BORA" ~ "03_BORA/",
                                                    analysis1 == "BORG" ~ "04_BORG/",
                                                    analysis1 == "CHIA" ~ "05_CHIA/",
                                                    analysis1 == "EASZ" ~ "06_EASZ/",
                                                    analysis1 == "KARA" ~ "07_KARA/",
                                                    analysis1 == "KETE" ~ "08_KETE/",
                                                    analysis1 == "MARC" ~ "09_MARC/",
                                                    analysis1 == "MARE" ~ "10_MARE/",
                                                    analysis1 == "NDAM" ~ "11_NDAM/",
                                                    analysis1 == "NGAN" ~ "12_NGAN/",
                                                    analysis1 == "ROMA" ~ "13_ROMA/",
                                                    analysis1 == "SHEK" ~ "14_SHEK/",
                                                    analysis1 == "SOMB" ~ "15_SOMB/",
                                                    .default = "16_combined/"),
                                   dplyr::case_when(stringr::str_detect(analysis1,
                                                                        "european") ~ "01_european-hybrids/",
                                                    stringr::str_detect(analysis1,
                                                                        "west") ~ "02_west-african-hybrids/",
                                                    stringr::str_detect(analysis1,
                                                                        "east") ~ "03_east-african-hybrids/",
                                                    .default = ""),
                                   software1,
                                   "-",
                                   data1,
                                   "-",
                                   stringr::str_replace_all(analysis1,
                                                            "_",
                                                            "-"),
                                   ".csv")) %>%
    dplyr::select(rs,
                  chr,
                  pos,
                  row_num,
                  mean_a,
                  mean_b,
                  mean_c,
                  z_a,
                  z_b,
                  z_c) %>%
    dplyr::rename(mean_a_1 = mean_a,
                  mean_b_1 = mean_b,
                  mean_c_1 = mean_c,
                  z_a_1 = z_a,
                  z_b_1 = z_b,
                  z_c_1 = z_c)
  # read input data 2
  input2 <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                   dplyr::case_when(software2 == "mosaic" ~ "01_mosaic/",
                                                    software2 == "elai" ~ "02_elai/"),
                                   "02_output/",
                                   dplyr::case_when(data2 == "hd" ~ "01_hd/",
                                                    data2 == "ld" ~ "02_ld/"),
                                   dplyr::case_when(analysis2 == "ALEN" ~ "01_ALEN/",
                                                    analysis2 == "ANKO" ~ "02_ANKO/",
                                                    analysis2 == "BORA" ~ "03_BORA/",
                                                    analysis2 == "BORG" ~ "04_BORG/",
                                                    analysis2 == "CHIA" ~ "05_CHIA/",
                                                    analysis2 == "EASZ" ~ "06_EASZ/",
                                                    analysis2 == "KARA" ~ "07_KARA/",
                                                    analysis2 == "KETE" ~ "08_KETE/",
                                                    analysis2 == "MARC" ~ "09_MARC/",
                                                    analysis2 == "MARE" ~ "10_MARE/",
                                                    analysis2 == "NDAM" ~ "11_NDAM/",
                                                    analysis2 == "NGAN" ~ "12_NGAN/",
                                                    analysis2 == "ROMA" ~ "13_ROMA/",
                                                    analysis2 == "SHEK" ~ "14_SHEK/",
                                                    analysis2 == "SOMB" ~ "15_SOMB/",
                                                    .default = "16_combined/"),
                                   dplyr::case_when(stringr::str_detect(analysis2,
                                                                        "european") ~ "01_european-hybrids/",
                                                    stringr::str_detect(analysis2,
                                                                        "west") ~ "02_west-african-hybrids/",
                                                    stringr::str_detect(analysis2,
                                                                        "east") ~ "03_east-african-hybrids/",
                                                    .default = ""),
                                   software2,
                                   "-",
                                   data2,
                                   "-",
                                   stringr::str_replace_all(analysis2,
                                                            "_",
                                                            "-"),
                                   ".csv")) %>%
    dplyr::select(rs,
                  chr,
                  pos,
                  row_num,
                  mean_a,
                  mean_b,
                  mean_c,
                  z_a,
                  z_b,
                  z_c) %>%
    dplyr::rename(mean_a_2 = mean_a,
                  mean_b_2 = mean_b,
                  mean_c_2 = mean_c,
                  z_a_2 = z_a,
                  z_b_2 = z_b,
                  z_c_2 = z_c)
  # join input data
  input <- dplyr::full_join(input1, input2)
  # plot ancestry a
  plot_a <- ggplot2::ggplot(input,
                            ggplot2::aes(x = mean_a_1,
                                         y = mean_a_2,
                                         colour = as.factor(chr),
                                         fill = as.factor(chr))) +
    ggplot2::geom_point(alpha = 0.05,
                        shape = 21,
                        stroke = 0.5) +
    ggplot2::geom_smooth(method = lm,
                         se = F) +
    viridis::scale_colour_viridis(direction = -1,
                                  discrete = T,
                                  name = "Chromosome") +
    viridis::scale_fill_viridis(direction = -1,
                                discrete = T,
                                name = "Chromosome") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6),
                    fill = ggplot2::guide_legend(ncol = 6)) +
    ggplot2::labs(x = paste0(stringr::str_to_upper(software1), " European *Bos taurus* ancestry"),
                  y = paste0(stringr::str_to_upper(software2), " European *Bos taurus* ancestry")) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                   axis.title.y = ggtext::element_markdown())
  # plot ancestry b
  plot_b <- ggplot2::ggplot(input,
                            ggplot2::aes(x = mean_b_1,
                                         y = mean_b_2,
                                         colour = as.factor(chr),
                                         fill = as.factor(chr))) +
    ggplot2::geom_point(alpha = 0.05,
                        shape = 21,
                        stroke = 0.5) +
    ggplot2::geom_smooth(method = lm,
                         se = F) +
    viridis::scale_colour_viridis(direction = -1,
                                  discrete = T,
                                  name = "Chromosome") +
    viridis::scale_fill_viridis(direction = -1,
                                discrete = T,
                                name = "Chromosome") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6),
                    fill = ggplot2::guide_legend(ncol = 6)) +
    ggplot2::labs(x = paste0(stringr::str_to_upper(software1), " African *Bos taurus* ancestry"),
                  y = paste0(stringr::str_to_upper(software2), " African *Bos taurus* ancestry")) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                   axis.title.y = ggtext::element_markdown())
  # plot ancestry c
  plot_c <- ggplot2::ggplot(input,
                            ggplot2::aes(x = mean_c_1,
                                         y = mean_c_2,
                                         colour = as.factor(chr),
                                         fill = as.factor(chr))) +
    ggplot2::geom_point(alpha = 0.05,
                        shape = 21,
                        stroke = 0.5) +
    ggplot2::geom_smooth(method = lm,
                         se = F) +
    viridis::scale_colour_viridis(direction = -1,
                                  discrete = T,
                                  name = "Chromosome") +
    viridis::scale_fill_viridis(direction = -1,
                                discrete = T,
                                name = "Chromosome") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6),
                    fill = ggplot2::guide_legend(ncol = 6)) +
    ggplot2::labs(x = paste0(stringr::str_to_upper(software1), " *Bos indicus* ancestry"),
                  y = paste0(stringr::str_to_upper(software2), " *Bos indicus* ancestry")) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                   axis.title.y = ggtext::element_markdown())
  # combine plots
  plot <- plot_a +
    patchwork::guide_area() +
    plot_b +
    plot_c +
    patchwork::plot_annotation(tag_levels = 'A') +
    patchwork::plot_layout(guides = "collect")
  # save plot
  png(paste0("07_local-ancestry-analysis/03_comparison/",
             software1,
             dplyr::case_when(software2 == software1 ~ "",
                              .default = paste0("-",
                                                software2)),
             "-",
             data1,
             dplyr::case_when(data2 == data1 ~ "",
                              .default = paste0("-",
                                                data2)),
             "-",
             stringr::str_replace_all(analysis1,
                                      "_",
                                      "-"),
             dplyr::case_when(analysis2 == analysis1 ~ "",
                              .default = paste0("-",
                                                stringr::str_replace_all(analysis2,
                                                                         "_",
                                                                         "-"))),
             "-correlation.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
}

### set analyses
analyses <- c("selected_european_hybrids",
              "selected_trypanotolerant_african_hybrids",
              "trypanosusceptible_african_hybrids")

### apply function to check correlation between mosaic and elai for hd and ld data
parallel::mclapply(analyses,
                   correlation_plot,
                   mc.cores = 15)
parallel::mclapply(analyses,
                   correlation_plot,
                   data1 = "ld",
                   mc.cores = 15)