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
#' @param facet Plot facets if multiple phenotypes and patterns provided (default = \code{"pheno"}).
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
                              facet = "pheno", ...) {
  plot_type <- match.arg(plot_type)
  
  x$patterns <- dplyr::filter(x$patterns,
                              max_lod >= min_lod)
  
  m <- x$patterns$founders %in% patterns
  patterns <- x$patterns$founders[m]
  pheno <- x$patterns$pheno[m]
  tmp <- x$scan
  colnames(tmp) <- pheno
  x$scan <- modify_object(x$scan, tmp[, m, drop = FALSE])
  
  pattern <- matrix(patterns, nrow(x$scan), ncol(x$scan), byrow = TRUE)
  
  tmp <- class(x$coef)
  x$coef <- x$coef[m]
  if(length(x$coef) > 1 & length(unique(pheno)) > 1) {
    names(x$coef) <- paste0(pheno, "_", names(x$coef))
  }
  class(x$coef) = tmp
  
  switch(plot_type,
         lod = autoplot(x$scan, map, lodcolumn = lodcolumn,
                        pattern = pattern, 
                        facet = facet, ...),
         coef = autoplot(x$coef, map, columns, ...),
         coef_and_lod = autoplot(x$coef, map, columns, 
                                 scan1_output = x$scan,
                                 lodcolumn = lodcolumn,
                                 pattern_lod = pattern, 
                                 facet_lod = facet, ...))
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

