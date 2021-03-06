#' Plot of object of class listof_scan1coeff
#'
#' Plot object of class \code{listof_scan1coeff}, which is a list of objects of class \code{scan1coef}.
#'
#' @param x object of class \code{listof_scan1coeff}
#'
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2]{insert_pseudomarkers}}.
#'
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
#' @param pattern Use phenotype names as pattern.
#'
#' @param ... arguments for \code{\link{ggplot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#' 
#' @return object of class \code{\link[ggplot2]{ggplot}}
#'
#' @export
#' @importFrom dplyr bind_rows
#' 
ggplot_listof_scan1coef <- function(x, map, columns = NULL, col = NULL,
                                  scan1_output = NULL,
                                  facet = "pattern",
                                  ...) {
  if(is.list(map))
    map1 <- map[[1]]
  else
    map1 <- map

  if(!is.null(scan1_output)) {
    # make scan1 lod rownames agree with first set in map
    rownames(scan1_output) <- names(map1)[seq_len(nrow(scan1_output))]
  }

  # Reform as one scan1coef object.
  
  # Bind together coef matrices.
  coefs <- dplyr::bind_rows(lapply(x, function(x) as.data.frame(unclass(x))),
                            .id = "pheno")
  pheno <- matrix(coefs$pheno, ncol = length(x))
  
  ## Put map names as rownames of coefs
  coefs <- as.matrix(coefs[,-1])
  rownames(coefs) <- c(outer(names(map1), names(x), paste, sep = "_"))
  map1 <- rep(map1, times = length(x))
  names(map1) <- rownames(coefs)
  if(is.list(map))
    map[[1]] <- map1
  else
    map <- map1
  
  # Get attributes and class right.
  attr(coefs, "sample_size") <- attr(x[[1]], "sample_size")
  attr(coefs, "SE") <- attr(x[[1]], "SE")
  class(coefs) <- class(x[[1]])

  ggplot_coef(coefs, map, columns, col, scan1_output,
            facet = facet,
            pattern = pheno,
            patterns = "all", ...)
}
#' @method autoplot listof_scan1coef
#' @export
#' @export autoplot.listof_scan1coef
#' @rdname ggplot_listof_scan1coef
#'
#' @importFrom ggplot2 autoplot
#'
autoplot.listof_scan1coef <- function(x, ...)
  ggplot_listof_scan1coef(x, ...)
