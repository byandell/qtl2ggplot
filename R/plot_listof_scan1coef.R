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
                                  xlim = NULL,
                                  ...) {
  # stretch out map
  map <- lapply(x, function(x) x$map)
  names(map) <- NULL
  map <- unlist(map)
  if(!is.null(scan1_output)) {
    # make scan1 lod rownames agree with first set in map
    rownames(scan1_output$lod) <- names(map)[seq_len(nrow(scan1_output$lod))]
  }

  # bind together coef matrices
  coefs <- dplyr::bind_rows(lapply(x, function(x) as.data.frame(x$coef)),
                            .id = "pheno")
  pheno <- matrix(coefs$pheno, ncol = length(x))
  coefs <- as.matrix(coefs[,-1])
  rownames(coefs) <- names(map)
  
  # Reform as one scan1coef object
  x <- x[[1]]
  x$map <- map
  x$coef <- coefs
  
  plot_coef(x, columns, col, scan1_output, 
            facet = facet, ylim=ylim_coef,
            pattern = pheno,
            patterns = "all", ...)
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
