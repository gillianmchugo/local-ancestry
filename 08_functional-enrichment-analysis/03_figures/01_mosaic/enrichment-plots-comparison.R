# enrichment plots comparison
## compare enrichment plots of local ancestry results
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to combine gprofiler enrichment plots
combine_gprofiler_plots <- function(software = "mosaic", data = "hd"){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # read, crop and annotate plots
  a <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "gprofiler-",
                                 software,
                                 "-",
                                 data,
                                 "-european-hybrids.png")) %>%
    magick::image_crop("2620x1540+80+120")
  a_label <- magick::image_append(c(magick::image_blank(100,
                                                        1540,
                                                        color = "white"),
                                    a)) %>%
    magick::image_annotate("European hybrid",
                           degrees = 270,
                           gravity = "west",
                           location = "+45+240",
                           size = 60) %>%
    magick::image_composite(magick::image_blank(5,
                                                1450,
                                                color = pal[15]),
                            gravity = "west",
                            offset = "+100+20") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[17]),
                            gravity = "east",
                            offset = "+110-506") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[9]),
                            gravity = "east",
                            offset = "+110+3") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[28]),
                            gravity = "east",
                            offset = "+110+512")
  b <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "gprofiler-",
                                 software,
                                 "-",
                                 data,
                                 "-trypanotolerant-african-hybrids.png")) %>%
    magick::image_crop("2620x1540+80+120")
  b_label <- magick::image_append(c(magick::image_blank(100,
                                                        1540,
                                                        color = "white"),
                                    b)) %>%
    magick::image_annotate("Trypanotolerant African hybrid",
                           degrees = 270,
                           gravity = "West",
                           location = "+45+420",
                           size = 60) %>%
    magick::image_composite(magick::image_blank(5,
                                                1450,
                                                color = pal[5]),
                            gravity = "west",
                            offset = "+100+20") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[17]),
                            gravity = "east",
                            offset = "+110-506") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[9]),
                            gravity = "east",
                            offset = "+110+3") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[28]),
                            gravity = "east",
                            offset = "+110+512")
  c <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "gprofiler-",
                                 software,
                                 "-",
                                 data,
                                 "-trypanosusceptible-african-hybrids.png")) %>%
    magick::image_crop("2620x1540+80+120")
  c_label <- magick::image_append(c(magick::image_blank(100,
                                                        1540,
                                                        color = "white"),
                                    c)) %>%
    magick::image_annotate("Trypanosusceptible African hybrid",
                           degrees = 270,
                           gravity = "West",
                           location = "+45+465",
                           size = 60) %>%
    magick::image_composite(magick::image_blank(5,
                                                1450,
                                                color = pal[26]),
                            gravity = "west",
                            offset = "+100+20") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[17]),
                            gravity = "east",
                            offset = "+110-506") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[9]),
                            gravity = "east",
                            offset = "+110+3") %>%
    magick::image_composite(magick::image_blank(5,
                                                486,
                                                color = pal[28]),
                            gravity = "east",
                            offset = "+110+512")
  # read, crop and add space to legend and axis titles
  legend <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                      dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                       software == "elai" ~ "02_elai/"),
                                      dplyr::case_when(data == "hd" ~ "01_hd/",
                                                       data == "ld" ~ "02_ld/"),
                                      "gprofiler-",
                                      software,
                                      "-",
                                      data,
                                      "-european-hybrids.png")) %>%
    magick::image_crop("2700x120")
  legend_space <- magick::image_append(c(magick::image_blank(20,
                                                             120,
                                                             color = "white"),
                                         legend))
  x <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "gprofiler-",
                                 software,
                                 "-",
                                 data,
                                 "-european-hybrids.png")) %>%
    magick::image_crop("2700x140+0+1660")
  x_space <- magick::image_append(c(magick::image_blank(20,
                                                        120,
                                                        color = "white"),
                                    x))
  y <- magick::image_read(paste0("08_functional-enrichment/03_figures/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "gprofiler-",
                                 software,
                                 "-",
                                 data,
                                 "-european-hybrids.png")) %>%
    magick::image_crop("80x1800")
  y_space <- magick::image_composite(magick::image_blank(80,
                                                         4880,
                                                         color = "white"),
                                     y,
                                     gravity = "west")
  # combine plots
  stack <- magick::image_append(c(legend_space,
                                  a_label,
                                  b_label,
                                  c_label,
                                  x_space),
                                stack = T)
  plot <- magick::image_append(c(y_space,
                                 stack))
  # save plot
  magick::image_write(plot,
                      paste0("08_functional-enrichment/03_figures/",
                             dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                              software == "elai" ~ "02_elai/"),
                             dplyr::case_when(data == "hd" ~ "01_hd/",
                                              data == "ld" ~ "02_ld/"),
                             "gprofiler-",
                             software,
                             "-",
                             data,
                             "-comparison.png"))
}

### apply function
combine_gprofiler_plots()
combine_gprofiler_plots(data = "ld")

combine_gprofiler_plots(software = "elai")
combine_gprofiler_plots(software = "elai",
                        data = "ld")