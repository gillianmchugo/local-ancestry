# structure plots labels
## add labels to plots of structure results
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

### read images and add label annotations
#### k = 2 - 3
magick::image_read("05_structure-analysis/03_figures/hd-combined-k-02-03.png") %>%
  magick::image_composite(magick::image_blank(545,
                                              10,
                                              color = pal[17]),
                          gravity = "southwest",
                          offset = "+275+230") %>%
  magick::image_composite(magick::image_blank(690,
                                              10,
                                              color = pal[15]),
                          gravity = "southwest",
                          offset = "+860+230") %>%
  magick::image_composite(magick::image_blank(470,
                                              10,
                                              color = pal[9]),
                          gravity = "southwest",
                          offset = "+1590+230") %>%
  magick::image_composite(magick::image_blank(990,
                                              10,
                                              color = pal[5]),
                          gravity = "south",
                          offset = "-100+230") %>%
  magick::image_composite(magick::image_blank(1615,
                                              10,
                                              color = pal[26]),
                          gravity = "southeast",
                          offset = "+645+230") %>%
  magick::image_composite(magick::image_blank(515,
                                              10,
                                              color = pal[28]),
                          gravity = "southeast",
                          offset = "+90+230") %>%
  magick::image_write("05_structure-analysis/03_figures/hd-combined-k-02-03-labels.png")

#### k = 2 - 9
magick::image_read("05_structure-analysis/03_figures/hd-combined-k-02-09.png") %>%
  magick::image_composite(magick::image_blank(555,
                                              10,
                                              color = pal[17]),
                          gravity = "southwest",
                          offset = "+230+230") %>%
  magick::image_composite(magick::image_blank(700,
                                              10,
                                              color = pal[15]),
                          gravity = "southwest",
                          offset = "+825+230") %>%
  magick::image_composite(magick::image_blank(480,
                                              10,
                                              color = pal[9]),
                          gravity = "southwest",
                          offset = "+1565+230") %>%
  magick::image_composite(magick::image_blank(1010,
                                              10,
                                              color = pal[5]),
                          gravity = "south",
                          offset = "-105+230") %>%
  magick::image_composite(magick::image_blank(1650,
                                              10,
                                              color = pal[26]),
                          gravity = "southeast",
                          offset = "+610+230") %>%
  magick::image_composite(magick::image_blank(525,
                                              10,
                                              color = pal[28]),
                          gravity = "southeast",
                          offset = "+45+230") %>%
  magick::image_write("05_structure-analysis/03_figures/hd-combined-k-02-09-labels.png")

#### k = 10 - 17
magick::image_read("05_structure-analysis/03_figures/hd-combined-k-10-17.png") %>%
  magick::image_composite(magick::image_blank(555,
                                              10,
                                              color = pal[17]),
                          gravity = "southwest",
                          offset = "+230+230") %>%
  magick::image_composite(magick::image_blank(700,
                                              10,
                                              color = pal[15]),
                          gravity = "southwest",
                          offset = "+825+230") %>%
  magick::image_composite(magick::image_blank(480,
                                              10,
                                              color = pal[9]),
                          gravity = "southwest",
                          offset = "+1565+230") %>%
  magick::image_composite(magick::image_blank(1010,
                                              10,
                                              color = pal[5]),
                          gravity = "south",
                          offset = "-105+230") %>%
  magick::image_composite(magick::image_blank(1650,
                                              10,
                                              color = pal[26]),
                          gravity = "southeast",
                          offset = "+610+230") %>%
  magick::image_composite(magick::image_blank(525,
                                              10,
                                              color = pal[28]),
                          gravity = "southeast",
                          offset = "+45+230") %>%
  magick::image_write("05_structure-analysis/03_figures/hd-combined-k-10-17-labels.png")