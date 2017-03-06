#' Plot QTL effects along chromosome
#'
#' Plot estimated QTL effects along a chromosomes.
#'
#' @param x Estimated QTL effects ("coefficients") as obtained from
#' \code{\link[qtl2scan]{scan1coef}}.
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
#' @param top_panel_prop If `scan1_output` provided, this gives the
#' proportion of the plot that is devoted to the top panel.
#'
#' @param center Center coefficients around 0 if \code{TRUE} (default)
#' 
#' @param CC use CC colors if \code{TRUE} (default if at least 8 columns of \code{coef} element of \code{x})
#' 
#' @param ylim_max max range for ylim (default is 0.1 and 99.9 percentile)
#' @param xlim x-axis limits. If \code{NULL}, we use the range of the plotted
#' coefficients.
#' 
#' @param maxpos,maxcol used for vertical line if maxpos is not \code{NULL} or \code{NA}
#' 
#' @param ... Additional graphics parameters.
#'
#' @export
#' @importFrom graphics layout par
#'
#' @details
#' \code{plot_coefCC()} is the same as \code{plot_coef()}, but forcing
#' \code{columns=1:8} and using the Collaborative Cross colors,
#' \code{\link[CCSanger]{CCcolors}}.
#'
#' @seealso \code{\link[CCSanger]{CCcolors}}, \code{\link{plot_scan1}}, \code{\link{plot_snpasso}}
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno[,1]
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # calculate coefficients for chromosome 7
#' library(qtl2scan)
#' coef <- scan1coef(probs[,7], pheno, addcovar=covar)
#'
#' # plot QTL effects
#' autoplot(coef, columns=1:3, col=c("slateblue", "violetred", "green3"))
plot_coef <-
    function(x, columns=NULL, col=NULL, scan1_output=NULL,
             gap=25, ylim=NULL,
             bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects", top_panel_prop=0.65, 
             center = TRUE, 
             CC = (ncol(x$coef) > 7),
             ylim_max = stats::quantile(x$coef, c(0.001, 0.999), na.rm = TRUE),
             xlim = NULL,
             maxpos = NULL, maxcol = 1,
             ...)
{
    if(!is.null(xlim)) {
      wh <- range(which(x$map >= xlim[1] & x$map <= xlim[2]))
      wh[1] <- max(1, wh[1] - 1)
      wh[2] <- min(length(x$map), wh[2] - 1)
      xlim <- x$map[wh]
    }
    if(!is.null(scan1_output)) { # call internal function for both coef and LOD
        return(plot_coef_and_lod(x, columns=columns, col=col, scan1_output=scan1_output,
                                 gap=gap, ylim=ylim, xlim=xlim,
                                 bgcolor=bgcolor, altbgcolor=altbgcolor,
                                 ylab="QTL effects", xaxt=NULL, top_panel_prop=top_panel_prop,
                                 center = center, ...))
    }

    # Set up CC colors if possible.
    all_columns <- NULL
    if(CC) {
      all_columns <- 1:8
      if(is.null(columns)) {
        columns <- 1:8
      }
      if(is.null(col)) {
        col <- CCSanger::CCcolors[columns]
      }
      if(is.null(names(col))) {
        names(col) <- dimnames(x$coef)[[2]][columns]
      } else {
        dimnames(x$coef)[[2]][columns] <- names(col)[seq_along(columns)]
      }
    }
    if(is.null(columns))
      columns <- 1:ncol(x$coef)
    if(is.null(all_columns))
      all_columns <- columns

    map <- x$map
    if(is.null(map)) stop("Input needs to contain a map")

    # Center coef on mean per locus if TRUE
    if(center) {
      col_mean <- apply(x$coef[, all_columns], 1, mean, na.rm=TRUE)
      x$coef[,columns] <- x$coef[,columns] - col_mean
    }

    if(is.null(ylim)) {
        ylim <- range(x$coef[,columns], na.rm=TRUE)
        d <- diff(ylim) * 0.02 # add 2% on either side
        ylim <- ylim + c(-d, d)
        ylim[1] <- max(min(ylim_max), ylim[1])
        ylim[2] <- min(max(ylim_max), ylim[2])
    }
    # Winsorize data limits
    tmp <- dimnames(x$coef)
    x$coef <- apply(x$coef, 2, 
                    function(x, ylim) pmin(pmax(x, ylim[1], na.rm = TRUE), 
                                           ylim[2], na.rm = TRUE),
                    ylim)
    
    names(x)[names(x)=="coef"] <- "lod" # switch coef -> lod for use with plot_scan1()

    plot_coef_internal <- function(x, columns, ylim, xlim, col, gap, 
                                   bgcolor, atlbgcolor, ylab, 
                                   legend.position = "right", 
                                   legend.title = "geno",
                                   maxpos = NULL,
                                   lodcolumn,
                                   ...) {
      lodcolumn <- columns # in case this is passed along from plot_coef_and_lod
      p <- plot_scan1(x, lodcolumn=lodcolumn, ylim=ylim, xlim=xlim,
                      col=col, gap=gap, 
                      bgcolor=bgcolor, altbgcolor=altbgcolor,
                      ylab=ylab, 
                      legend.position = legend.position,
                      legend.title = legend.title, 
                      ...)
      if(!is.null(maxpos)) {
        if(!is.na(maxpos))
          p <- p + ggplot2::geom_vline(xintercept = maxpos, 
                                       linetype=2,
                                       col = maxcol)
      }
      p
    }
    plot_coef_internal(x, columns, ylim, xlim, col, gap, 
                       bgcolor, atlbgcolor, ylab, 
                       maxpos = maxpos, maxcol = maxcol, ...)
}

#' @export
#' @rdname plot_coef
plot_coefCC <-
    function(x, scan1_output=NULL, gap=25, ylim=NULL,
             bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects", ...)
{
    dimnames(x$coef)[[2]][1:8] <- names(CCSanger::CCcolors)
    plot_coef(x, columns=1:8, col = CCSanger::CCcolors,
              scan1_output=scan1_output, gap=gap,
              ylim=ylim, bgcolor=bgcolor, altbgcolor=altbgcolor,
              ylab=ylab, ...)
}

#' @export
#' @export autoplot.scan1coef
#' @method autoplot scan1coef
#' @rdname plot_coef
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.scan1coef <-
    function(x, columns=NULL, col=NULL, scan1_output=NULL, gap=25, ylim=NULL,
             bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects", ...)
{
    plot_coef(x, columns=columns, col=col, scan1_output=scan1_output,
              gap=gap, ylim=ylim,
              bgcolor=bgcolor, altbgcolor=altbgcolor,
              ylab=ylab, ...)
}

#' @method plot scan1coef
#' @export
#' @export plot.scan1coef
#' @rdname plot_coef
#' 
plot.scan1coef <- function(x, ...)
  autoplot.scan1coef(x, ...)
