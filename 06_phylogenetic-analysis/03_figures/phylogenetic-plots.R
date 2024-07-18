# phylogenetic plots
## plot phylogenetic results from treemix
### install required packages
if (!require(ape)) install.packages("ape")
if (!require(devtools)) install.packages("devtools")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggh4x)) install.packages("ggh4x")
if (!require(ggnewscale)) install.packages("ggnewscale")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggtext)) install.packages("ggtext")
if (!require(ggtree)) devtools::install_github("YuLab-SMU/ggtree")
if (!require(khroma)) install.packages("khroma")
if (!require(tidytree)) install.packages("tidytree")
if (!require(treeio)) install.packages("treeio")
if (!require(viridis)) install.packages("viridis")

### functions to create plots
#### function to read tree data
read_tree <- function(m, data = "hd"){
  # read tree
  tree <- ggtree::read.tree(paste0("06_phylogenetic-analysis/02_output/02_bootstrap-replicates/",
                                   dplyr::case_when(data == "hd" ~ "01_hd/",
                                                    data == "ld" ~ "02_ld/"),
                                   stringr::str_pad(m,
                                                    2,
                                                    side = "left",
                                                    pad = "0"),
                                   "/",
                                   data,
                                   ".treeout.gz"))
  # read bootstraps
  bootstrap <- treeio::read.newick(paste0("06_phylogenetic-analysis/02_output/02_bootstrap-replicates/",
                                          dplyr::case_when(data == "hd" ~ "01_hd/",
                                                           data == "ld" ~ "02_ld/"),
                                          stringr::str_pad(m,
                                                           2,
                                                           pad = "0",
                                                           side = "left"),
                                          "/",
                                          data,
                                          "_outtree.newick"))
  # change bootstrap data format
  bootstrap_data <- treeio::as_tibble(bootstrap)
  # add bootstraps as branch lengths
  tree_data <- treeio::as_tibble(tree) %>%
    dplyr::mutate(label = dplyr::coalesce(label,
                                          as.character(bootstrap_data$branch.length)))
  # change tree data format
  treeio::as.phylo(tree_data)
}

#### function to draw tree
draw_tree <- function(tree, x2 = 0.3525, space = 0.03){
  # plot tree
  plot <- ggtree::ggtree(tree,
                         ladderize = F) +
    ggtree::geom_tippoint(ggplot2::aes(colour = label),
                          show.legend = F,
                          size = 3) +
    ggplot2::scale_colour_manual(values = c(pal[11], pal[18], pal[22], pal[26], pal[4],
                                            pal[14], pal[24], "#666666", pal[29], pal[17],
                                            pal[19], pal[25], pal[2], pal[8], pal[13],
                                            pal[12], pal[9], pal[7], pal[5], pal[30],
                                            pal[23], pal[15], pal[1], pal[3], pal[28])) +
    ggtree::geom_label2(ggplot2::aes(fill = as.numeric(label),
                                     label = label,
                                     subset = !isTip),
                        colour = "white",
                        hjust = 1.25,
                        label.padding = ggplot2::unit(0.1,
                                                      "lines"),
                        size = 3.25) +
    ggtree::geom_tiplab(colour = "grey30",
                        offset = 0.003) +
    ggplot2::labs(x = "Drift parameter") +
    ggtree::geom_strip(barsize = 0,
                       label = "atop(European, italic('Bos taurus'))",
                       offset = space,
                       offset.text = 0.003,
                       parse = T,
                       taxa1 = "HOLS",
                       taxa2 = "JERS") +
    ggtree::geom_strip(color = pal[17],
                       extend = 0.3,
                       offset = space,
                       taxa1 = "HOLS",
                       taxa2 = "JERS") +
    ggtree::geom_strip(barsize = 0,
                       label = "atop(European, hybrid)",
                       offset = space,
                       offset.text = 0.003,
                       parse = T,
                       taxa1 = "ALEN",
                       taxa2 = "MARE") +
    ggtree::geom_strip(color = pal[15],
                       extend = 0.3,
                       offset = space,
                       taxa1 = "ALEN",
                       taxa2 = "MARE") +
    ggtree::geom_strip(barsize = 0,
                       label = "italic('Bos indicus')",
                       offset = space,
                       offset.text = 0.003,
                       parse = T,
                       taxa1 = "THAR",
                       taxa2 = "NELO") +
    ggtree::geom_strip(color = pal[28],
                       extend = 0.3,
                       offset = space,
                       taxa1 = "THAR",
                       taxa2 = "NELO") +
    ggtree::geom_strip(barsize = 0,
                       label = "italic('Bos gaurus')",
                       offset = space,
                       offset.text = 0.003,
                       parse = T,
                       taxa1 = "GAUR",
                       taxa2 = "GAUR") +
    ggtree::geom_strip(color = "#666666",
                       extend = 0.3,
                       offset = space,
                       taxa1 = "GAUR",
                       taxa2 = "GAUR") +
    ggplot2::theme_light() +
    ggplot2::scale_fill_viridis_c(begin = 0,
                                  direction = -1,
                                  end = 0.85,
                                  limits = c(0,
                                             100),
                                  name = "Bootstrap\nvalue") +
    ggplot2::scale_x_continuous(expand = c(0,
                                           0),
                                limits = c(-0.0025,
                                           x2)) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(draw.ulim = F,
                                                    draw.llim = F)) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.box = "horizontal",
                   legend.position = c(0.05,
                                       0.84),
                   panel.grid = ggplot2::element_blank())
}

#### function to draw edges
draw_edges <- function(tree, edges, maxm = 0.5){
  tree +
    ggnewscale::new_scale_colour() +
    ggtree::geom_taxalink(data = edges,
                          mapping = ggplot2::aes(taxa1 = from,
                                                 taxa2 = to,
                                                 color = weight),
                          arrow = ggtree::arrow(length = ggplot2::unit(0.02,
                                                                       "npc"))) +
    ggplot2::scale_colour_viridis_c(begin = 0.5,
                                    direction = -1,
                                    end = 0.95,
                                    limits = c(0,
                                               maxm),
                                    name = "Migration\nweight",
                                    option = "inferno") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(draw.ulim = F,
                                                      draw.llim = F)) +
    ggplot2::theme(legend.position = c(0.1025,
                                       0.84))
}

#### function to save plot
save_tree <- function(plot, m, data = "hd"){
  # save plot
  png(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-tree.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-tree.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

#### function to calculate residuals adapted from treemix code
plot_resid_custom <- function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T, wcols = "r"){
  c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  m = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  names(c) = rownames(c)
  names(m) = rownames(m)
  o = read.table(pop_order, as.is = T, comment.char = "", quote = "")
  se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  mse = apply(se, 1, mean)
  mse = mean(mse)
  c = c[order(names(c)), order(names(c))]
  m = m[order(names(m)), order(names(m))]
  tmp = c -m
  toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
  for(i in 1:nrow(o)){
    for( j in 1:nrow(o)){
      if (o[i,1] %in% names(tmp) ==F){
        print(paste("not found", o[i,1]))
      }
      if (o[j,1] %in% names(tmp) ==F){
        print(paste("not found", o[j,1]))
      }
      toplot[i, j] = tmp[which(names(tmp)==o[i,1]), which(names(tmp)==o[j,1])]
    }
  }
  if (usemax){
    m1 = max(abs(toplot), na.rm = T)
    max = m1*1.02
    min = -(m1*1.02)
  }
  names(toplot) = o[,1]
  return(toplot)
}

#### function to draw residual plot
draw_resid <- function(m, data = "hd"){
  # set populations
  populations <- c("HOLS", "ANGU", "JERS",
                   "ROMA", "CHIA", "MARC", "MARE", "ALEN",
                   "MUTU", "LAGU", "NDAG",
                   "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                   "ANKO", "NGAN", "EASZ", "KARA", "BORA",
                   "THAR", "GIR", "NELO",
                   "GAUR")
  # set group labels
  group_labels <- c("European<br>*Bos taurus*",
                    "European<br>hybrid",
                    "African<br>*Bos taurus*",
                    "Trypanotolerant<br>African hybrid",
                    "Trypanosusceptible<br>African hybrid",
                    "*Bos<br>indicus*",
                    "*Bos<br>gaurus*")
  # set group levels
  group_levels <- c("European Bos taurus",
                    "European hybrid",
                    "African Bos taurus",
                    "Trypanotolerant African hybrid",
                    "Trypanosusceptible African hybrid",
                    "Bos indicus",
                    "Bos gaurus")
  # read residuals
  resid <- plot_resid_custom(paste0("06_phylogenetic-analysis/02_output/02_bootstrap-replicates/",
                                    dplyr::case_when(data == "hd" ~ "01_hd/",
                                                     data == "ld" ~ "02_ld/"),
                                    stringr::str_pad(m,
                                                     2,
                                                     pad = "0",
                                                     side = "left"),
                                    "/",
                                    data),
                             "06_phylogenetic-analysis/03_figures/poporder.txt")
  # set row names
  rownames(resid) <- populations
  # identify triangular matrix
  tri <- lower.tri(resid,
                   diag = T)
  # extract triangular matrix values
  resid_tri <- data.frame(col = colnames(resid)[col(resid)[tri]],
                          row = rownames(resid)[row(resid)[tri]],
                          value = resid[tri])
  # prepare data
  resid_tri <- resid_tri %>%
    dplyr::mutate(population1 = factor(col,
                                       levels = populations),
                  population2 = factor(row,
                                       levels = populations),
                  group1 = dplyr::case_when(population1 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                            population1 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                            population1 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                            population1 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                            population1 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                            population1 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus",
                                            population1 %in% c("GAUR") ~ "Bos gaurus"),
                  group2 = dplyr::case_when(population2 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                            population2 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                            population2 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                            population2 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                            population2 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                            population2 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus",
                                            population2 %in% c("GAUR") ~ "Bos gaurus"),
                  group1 = factor(group1,
                                  labels = group_labels,
                                  levels = group_levels),
                  group2 = factor(group2,
                                  labels = group_labels,
                                  levels = group_levels))
  ## set custom tick lengths
  ticks <- dplyr::tibble(population1 = populations,
                         length = c(rep(c(0.3,
                                          -0.3),
                                        12),
                                    0.3)) %>%
    dplyr::mutate(group1 = dplyr::case_when(population1 %in% c("HOLS", "ANGU", "JERS") ~ "European Bos taurus",
                                            population1 %in% c("ROMA", "CHIA", "MARC", "MARE", "ALEN") ~ "European hybrid",
                                            population1 %in% c("MUTU", "LAGU", "NDAG") ~ "African Bos taurus",
                                            population1 %in% c("NDAM", "BORG", "SOMB", "KETE", "SHEK") ~ "Trypanotolerant African hybrid",
                                            population1 %in% c("ANKO", "NGAN", "EASZ", "KARA", "BORA") ~ "Trypanosusceptible African hybrid",
                                            population1 %in% c("THAR", "GIR", "NELO") ~ "Bos indicus",
                                            population1 %in% c("GAUR") ~ "Bos gaurus"),
                  population2 = population1,
                  group2 = group1)
  # plot residuals
  plot <- ggplot2::ggplot(resid_tri,
                          ggplot2::aes(x = interaction(population1,
                                                       group1),
                                       y = interaction(population2,
                                                       group2),
                                       fill = value)) +
    ggplot2::geom_raster() +
    viridis::scale_fill_viridis(direction = -1,
                                limits = c(-max(abs(resid_tri$value)),
                                           max(abs(resid_tri$value))),
                                name = "Standard\nerror",
                                option = "cividis") +
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
    ggplot2::scale_x_discrete(expand = c(0,
                                         0),
                              guide = ggh4x::guide_axis_nested(n.dodge = 2)) +
    ggplot2::scale_y_discrete(expand = c(0,
                                         0),
                              guide = ggh4x::guide_axis_nested(),
                              limits = rev) +
    ggplot2::coord_equal(clip = "off",
                         expand = F,
                         ylim = c(0.5,
                                  25.5)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   ggh4x.axis.nestline.x = ggplot2::element_line(colour = c(pal[17],
                                                                            pal[15],
                                                                            pal[9],
                                                                            pal[5],
                                                                            pal[26],
                                                                            pal[28],
                                                                            "#666666"),
                                                                 linewidth = 0.5),
                   ggh4x.axis.nestline.y = ggplot2::element_line(colour = rev(c(pal[17],
                                                                                pal[15],
                                                                                pal[9],
                                                                                pal[5],
                                                                                pal[26],
                                                                                pal[28],
                                                                                "#666666")),
                                                                 linewidth = 0.5),
                   ggh4x.axis.nesttext.x = ggtext::element_markdown(colour = "black",
                                                                    size = ggplot2::rel(1)),
                   ggh4x.axis.nesttext.y = ggtext::element_markdown(colour = "black",
                                                                    size = ggplot2::rel(1),
                                                                    hjust = 0.5),
                   legend.background = ggplot2::element_blank(),
                   legend.position = c(0.925,
                                       0.84),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(0.8)),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  # save plot
  png(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-resid.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-resid.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  return plot
  return(plot)
}

#### function to combine residual and tree plots
combine_resid <- function(tree_plot, m, data = "hd", x2 = 0.3925){
  # draw residual plot
  resid_plot <- draw_resid(m,
                           data)
  # edit residual plot
  resid_plot <- resid_plot +
    ggplot2::coord_equal(clip = "on",
                         ylim = c(0.5,
                                  25.5)) +
    ggplot2::scale_x_discrete(expand = c(0,
                                         0),
                              guide = ggh4x::guide_axis_nested(angle = 90)) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(0.7)),
                   axis.ticks.x = ggplot2::element_line(),
                   ggh4x.axis.nesttext.x = ggplot2::element_blank(),
                   ggh4x.axis.nesttext.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(1)))
  # edit tree plot
  tree_plot_edit <- tree_plot +
    ggplot2::scale_x_continuous(expand = c(0,
                                           0),
                                limits = c(-0.0025,
                                           x2))
  # set layout
  layout <- "AAB
             AAC"
  # combine plots
  plot <- tree_plot_edit +
    patchwork::guide_area() +
    resid_plot +
    patchwork::plot_annotation(tag_levels = 'A') +
    patchwork::plot_layout(design = layout,
                           guides = "collect") &
    ggplot2::theme(legend.box = "horizontal")
  # save plot
  png(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-tree-resid.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("06_phylogenetic-analysis/03_figures/",
             data,
             "-m-",
             stringr::str_pad(m,
                              2,
                              side = "left",
                              pad = "0"),
             "-tree-resid.pdf"),
      height = 6,
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

### apply functions to data
#### hd 0
##### read tree
hd_0 <- read_tree(0)

##### rotate branches
hd_0_rotate <- ape::rotateConstr(hd_0,
                                 c("GAUR",
                                   "NELO", "GIR", "THAR",
                                   "BORA", "KARA", "EASZ", "NGAN", "SHEK",
                                   "ANKO", "KETE", "SOMB", "BORG", "NDAM",
                                   "NDAG", "LAGU", "MUTU",
                                   "MARE", "MARC", "CHIA", "ROMA", "ALEN",
                                   "JERS", "ANGU", "HOLS"))

##### draw tree
hd_0_tree <- draw_tree(hd_0_rotate) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "NDAM",
                     taxa2 = "KETE") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAM",
                     taxa2 = "KETE") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "SHEK",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "ANKO",
                     taxa2 = "ANKO") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "NGAN",
                     taxa2 = "BORA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NGAN",
                     taxa2 = "BORA")

##### save tree
save_tree(hd_0_tree,
          m = 0)

##### combine residual plot
combine_resid(hd_0_tree,
              m = 0)

#### hd 3
##### read tree
hd_3 <- read_tree(3)

##### rotate branches
hd_3_rotate <- ape::rotateConstr(hd_3,
                                 c("GAUR",
                                   "NELO", "GIR", "THAR",
                                   "BORA", "KETE", "KARA", "EASZ", "SHEK",
                                   "ANKO", "NGAN", "NDAM", "SOMB", "BORG",
                                   "NDAG", "LAGU", "MUTU",
                                   "MARE", "MARC", "CHIA", "ROMA", "ALEN",
                                   "JERS", "ANGU", "HOLS"))

##### draw tree
hd_3_rotate_tree <- draw_tree(hd_3_rotate) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "BORG",
                     taxa2 = "NDAM") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "BORG",
                     taxa2 = "NDAM") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "SHEK",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "KETE",
                     taxa2 = "KETE") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "ANKO",
                     taxa2 = "NGAN") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "EASZ",
                     taxa2 = "KARA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "EASZ",
                     taxa2 = "KARA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "BORA",
                     taxa2 = "BORA")

##### set migration edges
hd_3_edges <- dplyr::tibble(from = c(tidytree::MRCA(hd_3_rotate_tree,
                                                    "BORA",
                                                    "KARA"),
                                     tidytree::MRCA(hd_3_rotate_tree,
                                                    "LAGU",
                                                    "MUTU"),
                                     tidytree::nodeid(hd_3_rotate_tree,
                                                      "KARA")),
                            to = c(tidytree::nodeid(hd_3_rotate_tree,
                                                    "NGAN"),
                                   tidytree::nodeid(hd_3_rotate_tree,
                                                    "KETE"),
                                   tidytree::MRCA(hd_3_rotate_tree,
                                                  "SOMB",
                                                  "BORG")),
                            weight = c(0.460952,
                                       0.422615,
                                       0.478954))

##### draw migration edges
hd_3_rotate_tree_edges <- draw_edges(hd_3_rotate_tree,
                                     hd_3_edges)

##### save tree
save_tree(hd_3_rotate_tree_edges,
          m = 3)

##### combine residual plot
combine_resid(hd_3_rotate_tree_edges,
              m = 3)

#### hd 12
##### read tree
hd_12 <- read_tree(12)

##### rotate branches
hd_12_rotate <- ape::rotateConstr(hd_12,
                                  c("GAUR",
                                    "THAR", "GIR", "NELO",
                                    "BORA", "SHEK", "KETE", "ANKO", "KARA", "EASZ",
                                    "NGAN", "NDAM", "NDAG", "SOMB", "BORG",
                                    "LAGU", "MUTU",
                                    "MARE", "ROMA", "MARC", "CHIA", "ALEN",
                                    "JERS", "ANGU", "HOLS"))

##### draw tree
hd_12_rotate_tree <- draw_tree(hd_12_rotate,
                               x2 = 0.4825,
                               space = 0.041) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.041,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.041,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "SOMB",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "BORG",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "NDAG",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "NDAM",
                     taxa2 = "NDAM") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "NGAN",
                     taxa2 = "ANKO") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.041,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "NGAN",
                     taxa2 = "ANKO") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "KETE",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.041,
                     taxa1 = "BORA",
                     taxa2 = "BORA")

##### set migration edges
hd_12_edges <- dplyr::tibble(from = c(tidytree::MRCA(hd_12_rotate_tree,
                                                     "LAGU",
                                                     "MUTU"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "KARA",
                                                     "MUTU"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "KARA",
                                                     "NGAN"),
                                      tidytree::nodeid(hd_12_rotate_tree,
                                                       "ALEN"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "KARA",
                                                     "HOLS"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "NDAM",
                                                     "MUTU"),
                                      tidytree::nodeid(hd_12_rotate_tree,
                                                       "KETE"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "KARA",
                                                     "MUTU"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "SHEK",
                                                     "HOLS"),
                                      tidytree::MRCA(hd_12_rotate_tree,
                                                     "HOLS",
                                                     "MARE"),
                                      tidytree::nodeid(hd_12_rotate_tree,
                                                       "GIR"),
                                      tidytree::nodeid(hd_12_rotate_tree,
                                                       "ANKO")),
                             to = c(tidytree::nodeid(hd_12_rotate_tree,
                                                     "KETE"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "NGAN"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "ANKO"),
                                    tidytree::MRCA(hd_12_rotate_tree,
                                                   "NDAM",
                                                   "MUTU"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "NDAM"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "SHEK"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "BORG"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "BORA"),
                                    tidytree::MRCA(hd_12_rotate_tree,
                                                   "KARA",
                                                   "NGAN"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "GAUR"),
                                    tidytree::MRCA(hd_12_rotate_tree,
                                                   "CHIA",
                                                   "MARC"),
                                    tidytree::nodeid(hd_12_rotate_tree,
                                                     "SOMB")),
                             weight = c(0.494092,
                                        0.228541,
                                        0.492801,
                                        0.320188,
                                        0.170937,
                                        0.427069,
                                        0.380918,
                                        0.27374,
                                        0.706014,
                                        0.122396,
                                        0.0298064,
                                        0.436959))

##### draw migration edges
hd_12_rotate_tree_edges <- draw_edges(hd_12_rotate_tree,
                                      hd_12_edges,
                                      maxm = 0.75)

##### save tree
save_tree(hd_12_rotate_tree_edges,
          m = 12)

##### combine residual plot
combine_resid(hd_12_rotate_tree_edges,
              m = 12,
              x2 = 0.5425)

#### ld 0
##### read tree
ld_0 <- read_tree(0,
                  data = "ld")

##### rotate branches
ld_0_rotate <- ape::rotateConstr(ld_0,
                                 c("GAUR",
                                   "NELO", "GIR", "THAR",
                                   "BORA", "KARA", "EASZ", "NGAN", "ANKO",
                                   "SHEK", "KETE", "SOMB","BORG", "NDAM",
                                   "NDAG", "LAGU", "MUTU",
                                   "MARE", "MARC", "CHIA", "ROMA", "ALEN",
                                   "JERS", "ANGU", "HOLS"))

##### draw tree
ld_0_tree <- draw_tree(ld_0_rotate,
                       x2 = 0.3025) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "MUTU",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "NDAM",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAM",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "ANKO",
                     taxa2 = "BORA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "ANKO",
                     taxa2 = "BORA")

##### save tree
save_tree(ld_0_tree,
          m = 0,
          data = "ld")

##### combine residual plot
combine_resid(ld_0_tree,
              m = 0,
              x2 = 0.3425,
              data = "ld")

#### ld 3
##### read tree
ld_3 <- read_tree(3,
                  data = "ld")

##### rotate branches
ld_3_rotate <- ape::rotateConstr(ld_3,
                                 c("GAUR",
                                   "NELO", "GIR", "THAR",
                                   "BORA", "KARA", "KETE", "EASZ", "SOMB", "SHEK",
                                   "NGAN", "ANKO", "NDAM", "NDAG", "BORG",
                                   "LAGU", "MUTU",
                                   "MARE", "MARC", "CHIA", "ROMA", "ALEN",
                                   "JERS", "ANGU", "HOLS"))

##### draw tree
ld_3_rotate_tree <- draw_tree(ld_3_rotate,
                              x2 = 0.3025) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "BORG",
                     taxa2 = "BORG") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAG",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAM",
                     taxa2 = "NDAM") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "ANKO",
                     taxa2 = "NGAN") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "SHEK",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "SHEK",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "EASZ",
                     taxa2 = "EASZ") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "KETE",
                     taxa2 = "KETE") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "KARA",
                     taxa2 = "BORA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "KARA",
                     taxa2 = "BORA")

##### set migration edges
ld_3_edges <- dplyr::tibble(from = c(tidytree::MRCA(ld_3_rotate_tree,
                                                    "MUTU",
                                                    "LAGU"),
                                     tidytree::MRCA(ld_3_rotate_tree,
                                                    "NGAN",
                                                    "SHEK"),
                                     tidytree::MRCA(ld_3_rotate_tree,
                                                    "MUTU",
                                                    "LAGU")),
                            to = c(tidytree::nodeid(ld_3_rotate_tree,
                                                    "SOMB"),
                                   tidytree::nodeid(ld_3_rotate_tree,
                                                    "BORG"),
                                   tidytree::nodeid(ld_3_rotate_tree,
                                                    "KETE")),
                            weight = c(0.461115,
                                       0.499212,
                                       0.400278))

##### draw migration edges
ld_3_rotate_tree_edges <- draw_edges(ld_3_rotate_tree,
                                     ld_3_edges)

##### save tree
save_tree(ld_3_rotate_tree_edges,
          m = 3,
          data = "ld")

##### combine residual plot
combine_resid(ld_3_rotate_tree_edges,
              m = 3,
              x2 = 0.3425,
              data = "ld")

#### ld 11
##### read tree
ld_11 <- read_tree(11,
                   data = "ld")

##### rotate branches
ld_11_rotate <- ape::rotateConstr(ld_11,
                                  c("GAUR",
                                    "NELO", "GIR", "THAR",
                                    "KETE", "BORA", "KARA", "EASZ", "SHEK", "ANKO",
                                    "NGAN", "NDAM", "NDAG", "SOMB", "BORG",
                                    "LAGU", "MUTU",
                                    "MARE", "MARC", "CHIA", "ROMA", "ALEN",
                                    "JERS", "ANGU", "HOLS"))

##### draw tree
ld_11_rotate_tree <- draw_tree(ld_11_rotate,
                               x2 = 0.3225) +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(African, italic('Bos taurus'))",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "MUTU",
                     taxa2 = "LAGU") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanotolerant, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "SOMB",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "BORG",
                     taxa2 = "SOMB") +
  ggtree::geom_strip(color = pal[9],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAG",
                     taxa2 = "NDAG") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "NDAM",
                     taxa2 = "NDAM") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "ANKO",
                     taxa2 = "NGAN") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "SHEK",
                     taxa2 = "SHEK") +
  ggtree::geom_strip(barsize = 0,
                     label = "atop(Trypanosusceptible, 'African hybrid')",
                     offset = 0.03,
                     offset.text = 0.003,
                     parse = T,
                     taxa1 = "EASZ",
                     taxa2 = "BORA") +
  ggtree::geom_strip(color = pal[26],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "EASZ",
                     taxa2 = "BORA") +
  ggtree::geom_strip(color = pal[5],
                     extend = 0.3,
                     offset = 0.03,
                     taxa1 = "KETE",
                     taxa2 = "KETE")

##### set migration edges
ld_11_edges <- dplyr::tibble(from = c(tidytree::MRCA(ld_11_rotate_tree,
                                                     "LAGU",
                                                     "MUTU"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "KETE"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "KETE"),
                                      tidytree::MRCA(ld_11_rotate_tree,
                                                     "LAGU",
                                                     "NDAM"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "ALEN"),
                                      tidytree::MRCA(ld_11_rotate_tree,
                                                     "LAGU",
                                                     "NDAM"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "SOMB"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "SOMB"),
                                      tidytree::MRCA(ld_11_rotate_tree,
                                                     "LAGU",
                                                     "BORG"),
                                      tidytree::MRCA(ld_11_rotate_tree,
                                                     "LAGU",
                                                     "BORG"),
                                      tidytree::nodeid(ld_11_rotate_tree,
                                                       "GAUR")),
                             to = c(tidytree::nodeid(ld_11_rotate_tree,
                                                     "KETE"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "SOMB"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "BORG"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "BORA"),
                                    tidytree::MRCA(ld_11_rotate_tree,
                                                   "LAGU",
                                                   "NDAM"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "KARA"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "ANKO"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "SHEK"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "NGAN"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "EASZ"),
                                    tidytree::nodeid(ld_11_rotate_tree,
                                                     "NDAM")),
                             weight = c(0.491752,
                                        0.440039,
                                        0.385515,
                                        0.237458,
                                        0.342787,
                                        0.20709,
                                        0.356412,
                                        0.332856,
                                        0.244088,
                                        0.208092,
                                        0.111604))

##### draw migration edges
ld_11_rotate_tree_edges <- draw_edges(ld_11_rotate_tree,
                                      ld_11_edges)

##### save tree
save_tree(ld_11_rotate_tree_edges,
          m = 11,
          data = "ld")

##### combine residual plot
combine_resid(ld_11_rotate_tree_edges,
              m = 11,
              x2 = 0.3625,
              data = "ld")