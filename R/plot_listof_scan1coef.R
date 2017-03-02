#' Plot of object of class listof_scan1coeff
#'
#' Plot object of class \code{listof_scan1coeff}, which is a list of objects of class \code{scan1coef}.
#'
#' @param x object of class \code{listof_scan1coeff}
#' @param columns Vector of columns to plot
#'
#' @param col Vector of colors, same length as \code{columns}. If
#' NULL, some default choices are made.
#' 
#' @param scan1_output If provided, we make a two-panel plot with
#' coefficients on top and LOD scores below. Should have just one LOD
#' score column; if multiple, only the first is used.
#'
#' @param facet Plot facets if multiple phenotypes and group provided (default = \code{"pattern"}).
#' 
#' @param ylim_coef vertical limits for coef plot
#' 
#' @param ... arguments for \code{\link[qtl2scan]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @export
#' @importFrom dplyr bind_rows
plot_listof_scan1coef <- function(x, columns = NULL, col = NULL,
                                  scan1_output = NULL,
                                  facet = "pattern",
                                  ylim_coef = NULL,
                                  ...) {
  if(!is.null(scan1_output)) {
    return(plot_listof_coef_and_lod(x, columns, col, scan1_output, 
                                    facet, ylim_coef, ...))
  }
  
  # stretch out map
  map <- unlist(lapply(x, function(x) x$map))
  
  # bind together coef matrices
  coefs <- dplyr::bind_rows(lapply(x, function(x) as.data.frame(x$coef)),
                            .id = "pheno")
  pheno <- ordered(coefs$pheno, levels = names(x))
  coefs <- as.matrix(coefs[,-1])
  rownames(coefs) <- names(map)
  
  # Reform as one scan1coef object
  x <- x[[1]]
  x$map <- map
  x$coef <- coefs
  
  autoplot(x, columns, col, NULL, pattern = pheno,
       patterns = "all", facet = facet, ylim = ylim_coef, ...)
}
#' @method autoplot listof_scan1coef
#' @export
#' @export autoplot.listof_scan1coef
#' @rdname plot_listof_scan1coef
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.listof_scan1coef <- function(x, ...)
  plot_listof_scan1coef(x, ...)

#' @method plot listof_scan1coef
#' @export
#' @export plot.listof_scan1coef
#' @rdname plot_listof_scan1coef
#' 
plot.listof_scan1coef <- function(x, ...)
  autoplot.listof_scan1coef(x, ...)

#' @importFrom grid grid.layout grid.newpage pushViewport viewport
#' 
plot_listof_coef_and_lod <- function(x, columns, col, scan1_output, facet, ylim_coef,
                                     legend.position = "right", main = FALSE,
                                     maxpos = NULL, maxcol = 1,
                                     top_panel_prop = 0.65, ...) {
  
  if(is.null(maxpos)) { # include vertical line at max lod
    maxpos <- (dplyr::arrange(
      summary(scan1_output),
      dplyr::desc(lod)))$pos[1]
  }
  
  # 2 x 1 panels
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(nrow = 2,
                                 heights=c(top_panel_prop, 
                                           1-top_panel_prop))))
  
  p1 <- plot_listof_scan1coef(x, columns, col, NULL, facet, ylim_coef,
                              main = main, 
                              legend.position = legend.position,
                              maxpos = maxpos, maxcol = maxcol,
                              ...)
  print(p1,
        vp = grid::viewport(layout.pos.row = 1,
                            layout.pos.col = 1))
  
  p2 <- plot_scan1(scan1_output, 
                   legend.position = legend.position, ...)
  if(!is.na(maxpos))
    p2 <- p2 + ggplot2::geom_vline(xintercept = maxpos, 
                                   linetype=2,
                                   col = maxcol)
  print(p2,
        vp = grid::viewport(layout.pos.row = 2,
                            layout.pos.col = 1))
}
