#' Plot QTL effects along chromosome
#'
#' Plot estimated QTL effects along a chromosomes.
#'
#' @param x Estimated QTL effects ("coefficients") as obtained from
#' \code{\link[qtl2]{scan1coef}}.
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
#' @param gap Gap between chromosomes.
#'
#' @param ylim y-axis limits. If \code{NULL}, we use the range of the plotted
#' coefficients.
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param ylab y-axis label
#'
#' @param xlim x-axis limits. If \code{NULL}, we use the range of the plotted
#' coefficients.
#' 
#' @param colors Colors to use for plotting.
#'
#' @param ... Additional graphics parameters.
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @export
#' @importFrom graphics layout par
#' @importFrom grid grid.layout grid.newpage pushViewport viewport
#' @importFrom qtl2 align_scan1_map max_scan1
#' @importFrom dplyr arrange desc
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 theme element_blank geom_vline
#'
#' @details
#' \code{ggplot_coefCC()} is the same as \code{ggplot_coef()}, but forcing
#' \code{columns=1:8} and using the Collaborative Cross colors,
#' \code{\link[qtl2]{CCcolors}}.
#'
#' @seealso \code{\link{ggplot_scan1}}, \code{\link{ggplot_snpasso}}
#'
#' @examples
#' # load qtl2 package for data and genoprob calculation
#' library(qtl2)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' coef <- scan1coef(probs[,7], pheno, addcovar=covar)
#'
#' # plot QTL effects
#' library(ggplot2)
#' autoplot(coef, map[7], columns=1:3, col=c("slateblue", "violetred", "green3"))
ggplot_coef <-
    function(x, map, columns = NULL, col = NULL, scan1_output = NULL,
             gap = 25, ylim = NULL,
             bgcolor = "gray90", altbgcolor = "gray85",
             ylab = "QTL effects",
             xlim = NULL,
             ...)
{
    if(is.null(map)) stop("map is NULL")
    if(!is.list(map)) map <- list(" " = map) # if a vector, treat it as a list with no names
      
    # align scan1 output and map
    tmp <- qtl2::align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map
    
    if(nrow(x) != length(unlist(map)))
      stop("nrow(x) [", nrow(x), "] != number of positions in map [",
           length(unlist(map)), "]")
      
    if(!is.null(scan1_output)) { # call internal function for both coef and LOD
      return(ggplot_coef_and_lod(x, map, columns, col, scan1_output,
                                 gap, ylim, bgcolor, altbgcolor, ylab,
                                 xlim = xlim,
                                 ...))
    }

    ggplot_coef_internal(x, map, columns, ylim, xlim, col, gap,
                       bgcolor, altbgcolor = altbgcolor, ylab, ...)
}

ggplot_coef_internal <- function(x, map, columns, ylim, xlim, col, gap,
                               bgcolor, altbgcolor, ylab,
                               legend.position = "right",
                               legend.title = "geno",
                               maxpos = NULL, maxcol = 1,
                               lodcolumn,
                               scales = "free",
                               center = TRUE,
                               ylim_max = stats::quantile(x, c(0.001, 0.999), na.rm = TRUE),
                               colors = if(ncol(x) > 7) qtl2::CCcolors else NULL,
                               ...) {
  
  ## Replicate map as needed to be same size as x
  ## this is important if x came from plot_listof_scan1coef
  map <- rep(map, length = nrow(x))
  
  ## Change names and colors if used.
  all_columns <- NULL
  if(!is.null(colors)) {
    all_columns <- seq_along(colors)
    if(is.null(columns)) {
      columns <- seq_along(colors)
    }
    if(is.null(col)) {
      col <- colors[columns]
    }
    if(is.null(names(col))) {
      names(col) <- dimnames(x)[[2]][columns]
    } else {
      dimnames(x)[[2]][columns] <- names(col)[seq_along(columns)]
    }
  }
  if(is.null(columns))
    columns <- 1:ncol(x)
  if(is.null(all_columns))
    all_columns <- columns
  
  # Center coef on mean per locus if TRUE
  if(center) {
    col_mean <- apply(unclass(x)[, all_columns, drop=FALSE], 1, mean, na.rm=TRUE)
    x[, columns] <- unclass(x)[,columns,drop=FALSE] - col_mean
  }
  
  if(is.null(ylim)) {
    ylim <- range(unclass(x)[,columns], na.rm=TRUE)
    d <- diff(ylim) * 0.02 # add 2% on either side
    ylim <- ylim + c(-d, d)
    ylim[1] <- max(min(ylim_max), ylim[1])
    ylim[2] <- min(max(ylim_max), ylim[2])
  }
  # Winsorize data limits
  x <- apply(x, 2,
             function(x, ylim) pmin(pmax(x, ylim[1], na.rm = TRUE),
                                    ylim[2], na.rm = TRUE),
             ylim)
  
  if(is.list(map))
    map <- map[[1]]
  if(!is.null(xlim)) {
    wh <- range(which(map >= xlim[1] & map <= xlim[2]))
    wh[1] <- max(1, wh[1] - 1)
    wh[2] <- min(length(map), wh[2] - 1)
    xlim <- map[wh]
  }
  
  lodcolumn <- columns # in case this is passed along from plot_coef_and_lod
  p <- ggplot_scan1(x, map, lodcolumn=lodcolumn, ylim=ylim, xlim=xlim,
                  col=col, gap=gap,
                  bgcolor=bgcolor, altbgcolor = altbgcolor,
                  ylab=ylab,
                  legend.position = legend.position,
                  legend.title = legend.title,
                  scales = scales,
                  ...)
  if(!is.null(maxpos)) {
    if(!is.na(maxpos))
      p <- p + ggplot2::geom_vline(xintercept = maxpos,
                                   linetype=2,
                                   col = maxcol)
  }
  p
}

#' @export
#' @rdname ggplot_coef
ggplot_coefCC <- function(x, map, 
                          colors = qtl2::CCcolors, ...) {
  ncols <- seq_along(colors)
  assertthat::assert_that(dim(x)[2] >= length(ncols))
  
  dimnames(x)[[2]][ncols] <- names(colors)
  ggplot_coef(x, map, ncols, colors, ...)
}

#' @export
#' @export autoplot.scan1coef
#' @method autoplot scan1coef
#' @rdname ggplot_coef
#'
#' @importFrom ggplot2 autoplot
#'
autoplot.scan1coef <- function(x, ...) {
  ggplot_coef(x, ...)
}
