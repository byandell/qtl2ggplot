#' Plot scan pattern
#' 
#' @param x object of class \code{\link{scan_pattern}}
#' @param plot_type type of plot from \code{c("lod","coef")}
#' @param patterns allele patterns to plot (default all)
#' @param columns columns used for coef plot
#' @param min_lod minimum LOD peak for contrast to be retained
#' @param ylim_coef vertical limits for coef plot
#' @param ... additional parameters
#'
#' @export
#' 
#' @importFrom dplyr bind_cols filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes geom_path geom_vline ggplot ggtitle
#' 
plot_scan_pattern <- function(x,
                                  plot_type=c("lod","coef","coef_and_lod"),
                                  patterns = x$patterns$founders,
                                  columns = 1:3,
                                  min_lod = 3,
                                  ylim_coef = c(-2,2),
                                  ...) {
  plot_type <- match.arg(plot_type)
  
  x$patterns <- dplyr::filter(x$patterns,
                              max_lod >= min_lod)
  
  patterns <- patterns[patterns %in% x$patterns$founders]
  x$scan$lod <- x$scan$lod[,patterns, drop=FALSE]
  tmp <- class(x$coef)
  x$coef <- x$coef[patterns]
  class(x$coef) = tmp
  
  switch(plot_type,
         lod = plot(x$scan, lodcolumn = seq_along(patterns), ...),
         coef = plot(x$coef, columns, ylim = ylim_coef, ...),
         coef_and_lod = plot(x$coef, columns, scan1_output = x$scan))
}
#' @method autoplot scan_pattern
#' @export
#' @export autoplot.scan_pattern
#' @rdname plot_scan_pattern
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.scan_pattern <- function(x, ...)
  plot_scan_pattern(x, ...)
#' @export plot.scan_pattern
#' @export
#' @method plot scan_pattern
#' @rdname plot_scan_pattern
#' 
plot.scan_pattern <- function(x, ...) autoplot.scan_pattern(x, ...)
