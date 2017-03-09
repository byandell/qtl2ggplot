#' Plot scan pattern
#' 
#' @param x object of class \code{\link{scan_pattern}}
#' 
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2geno]{insert_pseudomarkers}}.
#'
#' @param plot_type type of plot from \code{c("lod","coef")}
#' @param patterns allele patterns to plot (default all)
#' @param columns columns used for coef plot
#' @param min_lod minimum LOD peak for contrast to be retained
#' @param lodcolumn columns used for scan1 plot (default all \code{patterns})
#' @param ... additional parameters
#'
#' @export
#' 
#' @importFrom dplyr bind_cols filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes geom_path geom_vline ggplot ggtitle
#' 
plot_scan_pattern <- function(x, map, plot_type = c("lod","coef","coef_and_lod"),
                              patterns = x$patterns$founders,
                              columns = 1:3,
                              min_lod = 3,
                              lodcolumn = seq_along(patterns),
                              ...) {
  plot_type <- match.arg(plot_type)
  
  x$patterns <- dplyr::filter(x$patterns,
                              max_lod >= min_lod)
  
  patterns <- patterns[patterns %in% x$patterns$founders]
  o <- order(-apply(x$scan[,patterns, drop=FALSE], 2, max))
  patterns <- patterns[o]
  
  sample_size <- attr(x$scan, "sample_size")
  x$scan <- x$scan[,patterns, drop=FALSE]
  attr(x$scan, "sample_size") <- sample_size
  class(x$scan) <- c("scan1", "matrix")

  tmp <- class(x$coef)
  x$coef <- x$coef[patterns]
  class(x$coef) = tmp
  
  switch(plot_type,
         lod = autoplot(x$scan, map, lodcolumn = lodcolumn, ...),
         coef = autoplot(x$coef, map, columns, ...),
         coef_and_lod = autoplot(x$coef, map, columns, 
                                 scan1_output = x$scan,
                                 lodcolumn = lodcolumn, ...))
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

#' @method plot scan_pattern
#' @export
#' @export plot.scan_pattern
#' @rdname plot_scan_pattern
#' 
plot.scan_pattern <- function(x, ...)
  autoplot.scan_pattern(x, ...)

