# load qtl2geno package for data and genoprob calculation
library(qtl2geno)

# read data
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

# insert pseudomarkers into map
map <- insert_pseudomarkers(iron$gmap, step=1)

# calculate genotype probabilities
probs <- calc_genoprob(iron, map, error_prob=0.002)

geno <- maxmarg(probs)

qtl2ggplot::plot_onegeno(geno, map)

qtl2ggplot::plot_onegeno(geno, map, 1:4)
CCcolors <- qtl2plot::CCcolors

source("R/plot_onegeno.R")

tmp <- get_geno_intervals(geno, map, 1:4)

# could also facet on chr and show individuals side by side
ggplot2::ggplot(dplyr::filter(tmp, is.na(ind))) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = chr, y = lo),
    col = "transparent") +
  ggplot2::ylab("Position (Mbp)") +
  ggplot2::xlab("Chromosome") +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = unclass(chr) - 0.25,
      xmax = unclass(chr) + 0.25, 
      ymin = X1,
      ymax = X2),
    fill = "white", col = "black") + 
  ggplot2::scale_fill_manual(values = CCcolors) +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = unclass(chr) - 0.25,
      xmax = unclass(chr) + 0.25, 
      ymin = lo,
      ymax = hi,
      fill = fill),
    dplyr::filter(tmp, !is.na(ind))) +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = unclass(chr) - 0.25,
      xmax = unclass(chr) + 0.25, 
      ymin = X1,
      ymax = X2),
    fill = "transparent", col = "black") +
  ggplot2::theme(legend.position = "none",
                 panel.background = ggplot2::element_rect(fill = "black"),
                 panel.grid.major.x = ggplot2::element_line(color = "grey80")) +
  ggplot2::scale_y_reverse() +
  ggplot2::facet_wrap(~ ind)

qtl2plot::plot_onegeno(geno, map)
