# local ancestry plots
## plot results of local ancestry analyses
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(khroma)) install.packages("khroma")
if (!require(parallel)) install.packages("parallel")
if (!require(patchwork)) install.packages("patchwork")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
if (!require(stringr)) install.packages("stringr")
if (!require(tibble)) install.packages("tibble")

### function to plot individual chromosomes
chromosome_plot <- function(chromosome, analysis, software = "mosaic", data = "hd"){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # plot results
  plot <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 "02_output/",
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                                                  analysis == "ANKO" ~ "02_ANKO/",
                                                  analysis == "BORA" ~ "03_BORA/",
                                                  analysis == "BORG" ~ "04_BORG/",
                                                  analysis == "CHIA" ~ "05_CHIA/",
                                                  analysis == "EASZ" ~ "06_EASZ/",
                                                  analysis == "KARA" ~ "07_KARA/",
                                                  analysis == "KETE" ~ "08_KETE/",
                                                  analysis == "MARC" ~ "09_MARC/",
                                                  analysis == "MARE" ~ "10_MARE/",
                                                  analysis == "NDAM" ~ "11_NDAM/",
                                                  analysis == "NGAN" ~ "12_NGAN/",
                                                  analysis == "ROMA" ~ "13_ROMA/",
                                                  analysis == "SHEK" ~ "14_SHEK/",
                                                  analysis == "SOMB" ~ "15_SOMB/",
                                                  .default = "16_combined/"),
                                 dplyr::case_when(stringr::str_detect(analysis,
                                                                      "european") ~ "01_european-hybrids/",
                                                  stringr::str_detect(analysis,
                                                                      "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                                                  stringr::str_detect(analysis,
                                                                      "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                                                  .default = ""),
                                 software,
                                 "-",
                                 data,
                                 "-",
                                 stringr::str_replace_all(analysis,
                                                          "_",
                                                          "-"),
                                 "-",
                                 stringr::str_pad(chromosome,
                                                  2,
                                                  pad = "0",
                                                  side = "left"),
                                 ".csv")) %>%
    ggplot2::ggplot(ggplot2::aes(x = pos)) +
    ggplot2::geom_area(ggplot2::aes(y = mean_a +
                                      mean_b +
                                      mean_c,
                                    fill = "a")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_b +
                                      mean_c,
                                    fill = "b")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_c,
                                    fill = "c")) +
    ggplot2::scale_fill_manual(labels = c(expression(paste("European ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste("African ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste(italic("Bos indicus"),
                                                           phantom("p")))),
                               values = c(pal[17],
                                          pal[9],
                                          pal[28])) +
    ggplot2::scale_x_continuous(expand=c(0,
                                         0),
                                labels = scales::label_number(scale_cut = c(0,
                                                                            "Mb" = 1000000)),
                                name = paste("Position on chromosome",
                                             chromosome)) +
    ggplot2::scale_y_continuous(expand=c(0,
                                         0),
                                name = "Mean ancestry") +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "top",
                   legend.text.align = 0,
                   legend.title = ggplot2::element_blank())
  # save plot
  png(paste0("07_local-ancestry-analysis/",
             dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                              software == "elai" ~ "02_elai/"),
             "03_figures/",
             dplyr::case_when(data == "hd" ~ "01_hd/",
                              data == "ld" ~ "02_ld/"),
             dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                              analysis == "ANKO" ~ "02_ANKO/",
                              analysis == "BORA" ~ "03_BORA/",
                              analysis == "BORG" ~ "04_BORG/",
                              analysis == "CHIA" ~ "05_CHIA/",
                              analysis == "EASZ" ~ "06_EASZ/",
                              analysis == "KARA" ~ "07_KARA/",
                              analysis == "KETE" ~ "08_KETE/",
                              analysis == "MARC" ~ "09_MARC/",
                              analysis == "MARE" ~ "10_MARE/",
                              analysis == "NDAM" ~ "11_NDAM/",
                              analysis == "NGAN" ~ "12_NGAN/",
                              analysis == "ROMA" ~ "13_ROMA/",
                              analysis == "SHEK" ~ "14_SHEK/",
                              analysis == "SOMB" ~ "15_SOMB/",
                              .default = "16_combined/"),
             dplyr::case_when(stringr::str_detect(analysis,
                                                  "european") ~ "01_european-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                              .default = ""),
             software,
             "-",
             data,
             "-",
             stringr::str_replace_all(analysis,
                                      "_",
                                      "-"),
             "-",
             stringr::str_pad(chromosome,
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
  pdf(paste0("07_local-ancestry-analysis/",
             dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                              software == "elai" ~ "02_elai/"),
             "03_figures/",
             dplyr::case_when(data == "hd" ~ "01_hd/",
                              data == "ld" ~ "02_ld/"),
             dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                              analysis == "ANKO" ~ "02_ANKO/",
                              analysis == "BORA" ~ "03_BORA/",
                              analysis == "BORG" ~ "04_BORG/",
                              analysis == "CHIA" ~ "05_CHIA/",
                              analysis == "EASZ" ~ "06_EASZ/",
                              analysis == "KARA" ~ "07_KARA/",
                              analysis == "KETE" ~ "08_KETE/",
                              analysis == "MARC" ~ "09_MARC/",
                              analysis == "MARE" ~ "10_MARE/",
                              analysis == "NDAM" ~ "11_NDAM/",
                              analysis == "NGAN" ~ "12_NGAN/",
                              analysis == "ROMA" ~ "13_ROMA/",
                              analysis == "SHEK" ~ "14_SHEK/",
                              analysis == "SOMB" ~ "15_SOMB/",
                              .default = "16_combined/"),
             dplyr::case_when(stringr::str_detect(analysis,
                                                  "european") ~ "01_european-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                              .default = ""),
             software,
             "-",
             data,
             "-",
             stringr::str_replace_all(analysis,
                                      "_",
                                      "-"),
             "-",
             stringr::str_pad(chromosome,
                              2,
                              pad = "0",
                              side = "left"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### function to plot round genome data
round_genome_plot <- function(analysis, size1 = 3, size2 = 2, software = "mosaic", data = "hd"){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # get input data
  input <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                  dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                   software == "elai" ~ "02_elai/"),
                                  "02_output/",
                                  dplyr::case_when(data == "hd" ~ "01_hd/",
                                                   data == "ld" ~ "02_ld/"),
                                  dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                                                   analysis == "ANKO" ~ "02_ANKO/",
                                                   analysis == "BORA" ~ "03_BORA/",
                                                   analysis == "BORG" ~ "04_BORG/",
                                                   analysis == "CHIA" ~ "05_CHIA/",
                                                   analysis == "EASZ" ~ "06_EASZ/",
                                                   analysis == "KARA" ~ "07_KARA/",
                                                   analysis == "KETE" ~ "08_KETE/",
                                                   analysis == "MARC" ~ "09_MARC/",
                                                   analysis == "MARE" ~ "10_MARE/",
                                                   analysis == "NDAM" ~ "11_NDAM/",
                                                   analysis == "NGAN" ~ "12_NGAN/",
                                                   analysis == "ROMA" ~ "13_ROMA/",
                                                   analysis == "SHEK" ~ "14_SHEK/",
                                                   analysis == "SOMB" ~ "15_SOMB/",
                                                   .default = "16_combined/"),
                                  dplyr::case_when(stringr::str_detect(analysis,
                                                                       "european") ~ "01_european-hybrids/",
                                                   stringr::str_detect(analysis,
                                                                       "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                                                   stringr::str_detect(analysis,
                                                                       "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                                                   .default = ""),
                                  software,
                                  "-",
                                  data,
                                  "-",
                                  stringr::str_replace_all(analysis,
                                                           "_",
                                                           "-"),
                                  ".csv"))
  # get chromosome lengths
  length <- input %>%
    dplyr::count(chr)
  # set label positions
  label_pos <- c()
  for(i in 1:29){
    label_pos <- c(label_pos,
                   sum(length$n[0:(i-1)],
                       length$n[i]/2))
  }
  # plot results
  plot <- ggplot2::ggplot(input,
                          ggplot2::aes(x = row_num)) +
    ggplot2::geom_area(ggplot2::aes(y = mean_a +
                                      mean_b +
                                      mean_c,
                                    fill = "a")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_b +
                                      mean_c,
                                    fill = "b")) +
    ggplot2::geom_area(ggplot2::aes(y = mean_c,
                                    fill = "c")) +
    ggplot2::annotate("text",
                      x = label_pos,
                      y = -0.1,
                      label = c(1:29),
                      size = c(rep(size1,
                                   20),
                               rep(size2,
                                   9)),
                      colour = "grey30") +
    ggplot2::scale_fill_manual(labels = c(expression(paste("European ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste("African ",
                                                           italic("Bos taurus"),
                                                           phantom("p"))),
                                          expression(paste(italic("Bos indicus"),
                                                           phantom("p")))),
                               values = c(pal[17],
                                          pal[9],
                                          pal[28])) +
    ggplot2::scale_x_continuous(expand=c(0,
                                         0)) +
    ggplot2::scale_y_continuous(expand=c(0,
                                         0),
                                limits = c(-1,
                                           1.01)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.text.align = 0,
                   legend.title = ggplot2::element_blank()) +
    ggplot2::coord_polar()
  # add line segments between chromosomes
  for(i in 1:29){
    plot <- plot +
      ggplot2::geom_segment(x = sum(length$n[1:i],
                                    1),
                            xend = sum(length$n[1:i],
                                       1),
                            y = 0,
                            yend = 1,
                            colour = "grey70",
                            linewidth = ggplot2::rel(0.25))
  }
  # return plot
  plot
}

### function to save single round genome plot
save_round_genome_plot <- function(analysis, software = "mosaic", data = "hd"){
  # get plot
  plot <- get(paste0(data,
                     "_",
                     analysis,
                     "_round_plot"))
  # save plot
  png(paste0("07_local-ancestry-analysis/",
             dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                              software == "elai" ~ "02_elai/"),
             "03_figures/",
             dplyr::case_when(data == "hd" ~ "01_hd/",
                              data == "ld" ~ "02_ld/"),
             dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                              analysis == "ANKO" ~ "02_ANKO/",
                              analysis == "BORA" ~ "03_BORA/",
                              analysis == "BORG" ~ "04_BORG/",
                              analysis == "CHIA" ~ "05_CHIA/",
                              analysis == "EASZ" ~ "06_EASZ/",
                              analysis == "KARA" ~ "07_KARA/",
                              analysis == "KETE" ~ "08_KETE/",
                              analysis == "MARC" ~ "09_MARC/",
                              analysis == "MARE" ~ "10_MARE/",
                              analysis == "NDAM" ~ "11_NDAM/",
                              analysis == "NGAN" ~ "12_NGAN/",
                              analysis == "ROMA" ~ "13_ROMA/",
                              analysis == "SHEK" ~ "14_SHEK/",
                              analysis == "SOMB" ~ "15_SOMB/",
                              .default = "16_combined/"),
             dplyr::case_when(stringr::str_detect(analysis,
                                                  "european") ~ "01_european-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                              .default = ""),
             software,
             "-",
             data,
             "-",
             stringr::str_replace_all(analysis,
                                      "_",
                                      "-"),
             "-round.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("07_local-ancestry-analysis/",
             dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                              software == "elai" ~ "02_elai/"),
             "03_figures/",
             dplyr::case_when(data == "hd" ~ "01_hd/",
                              data == "ld" ~ "02_ld/"),
             dplyr::case_when(analysis == "ALEN" ~ "01_ALEN/",
                              analysis == "ANKO" ~ "02_ANKO/",
                              analysis == "BORA" ~ "03_BORA/",
                              analysis == "BORG" ~ "04_BORG/",
                              analysis == "CHIA" ~ "05_CHIA/",
                              analysis == "EASZ" ~ "06_EASZ/",
                              analysis == "KARA" ~ "07_KARA/",
                              analysis == "KETE" ~ "08_KETE/",
                              analysis == "MARC" ~ "09_MARC/",
                              analysis == "MARE" ~ "10_MARE/",
                              analysis == "NDAM" ~ "11_NDAM/",
                              analysis == "NGAN" ~ "12_NGAN/",
                              analysis == "ROMA" ~ "13_ROMA/",
                              analysis == "SHEK" ~ "14_SHEK/",
                              analysis == "SOMB" ~ "15_SOMB/",
                              .default = "16_combined/"),
             dplyr::case_when(stringr::str_detect(analysis,
                                                  "european") ~ "01_european-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanotolerant") ~ "02_trypanotolerant-african-hybrids/",
                              stringr::str_detect(analysis,
                                                  "trypanosusceptible") ~ "03_trypanosusceptible-african-hybrids/",
                              .default = ""),
             software,
             "-",
             data,
             "-",
             stringr::str_replace_all(analysis,
                                      "_",
                                      "-"),
             "-round.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

### apply functions to data
#### set list of analyses
individual_plots <- c("ROMA", "CHIA", "MARC", "MARE", "ALEN",
                      "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                      "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                      "european_hybrids",
                      "trypanotolerant_african_hybrids",
                      "trypanosusceptible_african_hybrids",
                      "selected_european_hybrids",
                      "selected_trypanotolerant_african_hybrids")

#### save chromosome plots for all hybrids and comparisons
for (i in individual_plots){
  parallel::mclapply(1:29,
                     chromosome_plot,
                     analysis = i,
                     mc.cores = 15)
  parallel::mclapply(1:29,
                     chromosome_plot,
                     analysis = i,
                     data = "ld",
                     mc.cores = 15)
}

#### draw individual round genome plots for all hybrids and comparisons
for (i in individual_plots){
  assign(paste0("hd_",
                i,
                "_round_plot"),
         round_genome_plot(i))
  assign(paste0("ld_",
                i,
                "_round_plot"),
         round_genome_plot(i,
                           data = "ld"))
}

#### save individual round genome plots for all hybrids and comparisons
parallel::mclapply(individual_plots,
                   save_round_genome_plot,
                   mc.cores = 15)
parallel::mclapply(individual_plots,
                   save_round_genome_plot,
                   data = "ld",
                   mc.cores = 15)