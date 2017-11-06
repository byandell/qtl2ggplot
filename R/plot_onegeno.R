#' Plot one individual's genome-wide genotypes
#'
#' Plot one individual's genome-wide genotypes
#'
#' @param geno Imputed phase-known genotypes, as a list of matrices
#'     (as produced by \code{\link[qtl2geno]{maxmarg}}) or a list of
#'     three-dimensional arrays (as produced by \code{\link[qtl2geno]{guess_phase}}).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosomes to plot; a vector of character strings.
#' @param col Vector of colors for the different genotypes.
#' @param na_col Color for missing segments.
#' @param border Color of outer border around chromosome rectangles.
#' @param shift If TRUE, shift the chromosomes so they all start at 0.
#' @param bgcolor Color for background rectangle
#' @param chrwidth Total width of rectangles for each chromosome, as a
#'     fraction of the distance between them.
#' @param ... Additional graphics parameters
#'
#' @export
#' @importFrom graphics plot rect par axis title abline box
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # inferred genotype at a 28.6 cM on chr 16
#' geno <- maxmarg(probs)
#'
#' plot_onegeno(geno, map, shift = TRUE)
#' 
#' plot_onegeno(geno, map, ind=1:4)
plot_onegeno <-
    function(geno, map, ind=1, chr=NULL, col=NULL, na_col="white",
             border="black", shift=FALSE, bgcolor="gray90",
             chrwidth=0.5, ...)
{
    if(is.null(map)) stop("map is NULL")

    stopifnot(chrwidth > 0 && chrwidth < 1)

    # find common chr
    chr_geno <- names(geno)
    chr_map <- names(map)
    common_chr <- chr_geno[chr_geno %in% chr_map]
    if(length(common_chr) == 0)
        stop("No chr in common between geno and map")
    geno <- geno[common_chr]
    map <- map[common_chr]

    # subset chr if necessary
    if(!is.null(chr)) {
        if(any(!(chr %in% common_chr))) {
            chr <- chr[chr %in% common_chr]
            if(length(chr) == 0)
                stop("No chromosomes in common between geno, map, and chr")
            warning("Dropping some chr not found in geno and/or map")
        }
        geno <- geno[chr]
        map <- map[chr]
    }

    # same numbers of markers?
    nmar_map <- sapply(map, length)
    nmar_geno <- sapply(geno, ncol)
    if(any(nmar_geno != nmar_map))
        stop("Mismatch between numbers of markers between geno and map on chr ",
             paste(names(geno)[nmar_geno != nmar_map], collapse=", "))
    for(i in seq_along(geno)) {
        if(!all(names(map[[i]]) == colnames(geno[[i]])))
            stop("Mismatch between marker names on chr ", names(geno)[i])
    }

    # shift map to start at 0
    if(shift) map <- lapply(map, function(a) a-min(a,na.rm=TRUE))
    
    plot_onegeno_internal <-
        function(geno, map, ind, col=NULL, na_col="white",
                 border="black", bgcolor="gray90",
                 chrwidth=0.5,
                 xlab="Chromosome", ylab="Position (Mbp)",
                 ylim=NULL, main="",
                 vlines.col="gray80",
                 ...)
    {
        dots <- list(...)
        if(is.null(ylim)) ylim <- rev(range(unlist(map), na.rm=TRUE))

        nchr <- length(map)

        intervals <- get_geno_intervals(geno, map, ind)
        
        dims <- dplyr::mutate(
          dplyr::bind_rows(
            purrr::map(map, function(x) data.frame(t(range(x)))),
            .id = "chr"),
          chr = factor(chr, chr))
        
        chrwidth <- chrwidth / 2
        
        ## initial plot setup
        p <- ggplot2::ggplot(dims) +
          ggplot2::geom_point(
            ggplot2::aes(
              x = chr, y = X1),
            col = "transparent") +
          ggplot2::ylab(xlab) +
          ggplot2::xlab(ylab) +
          ggplot2::geom_rect(
            ggplot2::aes(
              xmin = unclass(chr) - chrwidth,
              xmax = unclass(chr) + chrwidth, 
              ymin = X1,
              ymax = X2),
            fill = na_col, col = border) +
          ggplot2::theme(
            legend.position = "none",
            panel.background = ggplot2::element_rect(fill = bgcolor))
        
        p <- p +
          ggplot_grid_lines(p, vlines.col=vlines.col, ...)

        # grid lines
        p <- p + ggplot_grid_lines(p, vlines, hlines, onechr)
        # grid lines
        if(!(length(vlines)==1 && is.na(vlines))) {
            if(is.null(vlines)) vlines <- 1:nchr
            abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
        }
        if(!(length(hlines)==1 && is.na(hlines))) {
            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
        }

        # x and y axis labels
        title(xlab=xlab, mgp=mgp.x)
        title(ylab=ylab, mgp=mgp.y)

        max_geno <- max(unlist(geno), na.rm=TRUE)
        if(is.null(col)) {
            if(max_geno <= 8) {
                col <- qtl2plot::CCcolors
            }
            else {
                warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
                col <- rep(qtl2plot::CCcolors, max_geno)
            }
        }
        else if(max_geno > length(col)) {
            warning("not enough colors; recycling them")
            col <- rep(col, max_geno)
        }

        for(i in 1:nchr) {
            g <- geno[[i]]

            # if completely missing the second chr, treat as if we have just the one
            this_chrwidth <- chrwidth
            if(!is.matrix(g) && all(is.na(g[,,2]))) {
                g <- rbind(g[,,1]) # make it a row matrix
                this_chrwidth <- this_chrwidth/2
            }

            if(is.matrix(g)) { # phase-known
                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i+this_chrwidth/2, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)

                addgenorect(g[1,], map[[i]], i-this_chrwidth/2, i+this_chrwidth/2,
                            col=col)

                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i+this_chrwidth/2, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
            }
            else {
                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)

                addgenorect(g[1,,1], map[[i]], i-this_chrwidth/2, i,
                            col=col)
                addgenorect(g[1,,2], map[[i]], i+this_chrwidth/2, i,
                            col=col)

                rect(i-this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map[[i]], na.rm=TRUE),
                     i, max(map[[i]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
            }
        }

        box()
    }

    plot_onegeno_internal(geno, map, ind, 
                          col=col, na_col=na_col, border=border,
                          bgcolor=bgcolor, chrwidth=chrwidth, ...)


}

get_geno_intervals <- function(geno, map, ind = 1) {
  # set up genotype intervals
  # for now, geno is reduced in parent function to one interval
  # want to be able to do this for multiple individuals.
  if(is.matrix(geno[[1]])){
    tmpfn <- function(x, ind) 
      geno2intervals(x$geno[ind,, drop = FALSE], x$map)
  } else {
    if(!is.array(geno[[i]]) || length(dim(geno[[i]])) != 3 ||
       dim(geno[[i]])[3] != 2)
      stop("geno should be an array individuals x positions x 2 haplotypes")

    tmpfn <- function(x, ind) 
      list(
        geno2intervals(x$geno[,,1][ind,, drop = FALSE], x$map),
        geno2intervals(x$geno[,,2][ind,, drop = FALSE], x$map))
  }

  dplyr::mutate(
    dplyr::bind_rows(
      purrr::map(
        purrr::transpose(
          purrr::map(
            purrr::transpose(
              list(geno = geno,
                   map = map)),
            tmpfn,
            ind)),
        dplyr::bind_rows,
        .id = "chr"),
      .id = "ind"),
    chr = factor(chr, names(map)[names(map) %in% chr]),
    col = factor(names(CCcolors)[geno], names(CCcolors)))
}

# convert vector of integer genotypes to intervals with common genotypes
# (start, end, genotype)
geno2intervals <- function(geno, map) {
  apply(geno, 1, geno2intervals_one, map)
}
geno2intervals_one <-
    function(geno, map)
{
    if(all(is.na(geno))) return(NULL)

    stopifnot(length(geno) == length(map))

    # drop missing values
    map <- map[!is.na(geno)]
    geno <- geno[!is.na(geno)]

    d <- diff(geno)
    xo_int <- which(d != 0)

    data.frame(lo=map[c(1,xo_int+1)],
               hi=map[c(xo_int, length(map))],
               geno=geno[c(xo_int, length(map))])

}
