# combine local ancestry plots
## combine local ancestry plots to compare analyses
### install required packages
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")

### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

### generate colour palette
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

### read legend
legend <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/01_hd/16_combined/01_european-hybrids/mosaic-hd-european-hybrids-01.png") %>%
  magick::image_crop("2700x170")

### read hd mosaic local ancestry plots
mosaic_hd_a <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/01_hd/16_combined/01_european-hybrids/mosaic-hd-european-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("HD",
                         gravity = "west",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              1400,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0")

mosaic_hd_b <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/01_hd/16_combined/02_trypanotolerant-african-hybrids/mosaic-hd-trypanotolerant-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

mosaic_hd_c <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/01_hd/16_combined/03_trypanosusceptible-african-hybrids/mosaic-hd-trypanosusceptible-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

### combine hd mosaic local ancestry plots
mosaic_hd <- magick::image_append(c(mosaic_hd_a,
                                    mosaic_hd_b,
                                    mosaic_hd_c))

### read ld mosaic local ancestry plots
mosaic_ld_a <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/02_ld/16_combined/01_european-hybrids/mosaic-ld-european-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("LD",
                         gravity = "west",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              1400,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0")

mosaic_ld_b <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/02_ld/16_combined/02_trypanotolerant-african-hybrids/mosaic-ld-trypanotolerant-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

mosaic_ld_c <- magick::image_read("07_local-ancestry-analysis/01_mosaic/03_figures/02_ld/16_combined/03_trypanosusceptible-african-hybrids/mosaic-ld-trypanosusceptible-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

### combine ld mosaic local ancestry plots
mosaic_ld <- magick::image_append(c(mosaic_ld_a,
                                    mosaic_ld_b,
                                    mosaic_ld_c))

### combine hd and ld mosaic local ancestry plots
mosaic <- magick::image_append(c(mosaic_hd,
                                 mosaic_ld),
                               stack = T)

### add mosaic label
mosaic_label <- magick::image_append(c(magick::image_blank(200,
                                                           3600,
                                                           color = "white"),
                                       mosaic)) %>%
  magick::image_annotate("MOSAIC",
                         degrees = 270,
                         gravity = "west",
                         location = "+75+205",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              3000,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0") %>%
  magick::image_scale("2700")

### read hd elai local ancestry plots
elai_hd_a <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/01_hd/16_combined/01_european-hybrids/elai-hd-european-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("HD",
                         gravity = "west",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              1400,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0")

elai_hd_b <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/01_hd/16_combined/02_trypanotolerant-african-hybrids/elai-hd-trypanotolerant-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

elai_hd_c <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/01_hd/16_combined/03_trypanosusceptible-african-hybrids/elai-hd-trypanosusceptible-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200")

### combine hd elai local ancestry plots
elai_hd <- magick::image_append(c(elai_hd_a,
                                  elai_hd_b,
                                  elai_hd_c))

### read ld elai local ancestry plots and add labels
elai_ld_a <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/02_ld/16_combined/01_european-hybrids/elai-ld-european-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("LD",
                         gravity = "west",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              1400,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0") %>%
  magick::image_annotate("European hybrid",
                         gravity = "south",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(1400,
                                              10,
                                              color = pal[[15]]),
                          gravity = "south",
                          offset = "+0+125")

elai_ld_b <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/02_ld/16_combined/02_trypanotolerant-african-hybrids/elai-ld-trypanotolerant-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("Trypanotolerant African hybrid",
                         gravity = "south",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(1400,
                                              10,
                                              color = pal[[5]]),
                          gravity = "south",
                          offset = "+0+125")

elai_ld_c <- magick::image_read("07_local-ancestry-analysis/02_elai/03_figures/02_ld/16_combined/03_trypanosusceptible-african-hybrids/elai-ld-trypanosusceptible-african-hybrids-round.png") %>%
  magick::image_crop("1800x1800+200") %>%
  magick::image_annotate("Trypanosusceptible African hybrid",
                         gravity = "south",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(1400,
                                              10,
                                              color = pal[[26]]),
                          gravity = "south",
                          offset = "+0+125")

### combine ld elai local ancestry plots
elai_ld <- magick::image_append(c(elai_ld_a,
                                  elai_ld_b,
                                  elai_ld_c))

### combine hd and ld elai local ancestry plots
elai <- magick::image_append(c(elai_hd,
                               elai_ld),
                             stack = T)

### add elai label
elai_label <- magick::image_append(c(magick::image_blank(200,
                                                         3600,
                                                         color = "white"),
                                     elai)) %>%
  magick::image_annotate("ELAI",
                         degrees = 270,
                         gravity = "west",
                         location = "+75+100",
                         size = 100) %>%
  magick::image_composite(magick::image_blank(10,
                                              3000,
                                              color = "grey70"),
                          gravity = "west",
                          offset = "+150+0") %>%
  magick::image_scale("2700")

### combine all plots and legend
mosaic_elai <- magick::image_append(c(legend,
                                      mosaic_label,
                                      elai_label),
                                    stack = T)

### save plot
magick::image_write(mosaic_elai,
                    "07_local-ancestry-analysis/03_comparison/local-ancestry-plots-comparison.png")