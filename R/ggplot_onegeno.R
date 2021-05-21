#' Plot one individual's genome-wide genotypes
#'
#' Plot one individual's genome-wide genotypes
#'
#' @param geno Imputed phase-known genotypes, as a list of matrices
#'     (as produced by \code{\link[qtl2]{maxmarg}}) or a list of
#'     three-dimensional arrays (as produced by \code{\link[qtl2]{guess_phase}}).
#' @param map Marker map (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosomes to plot; a vector of character strings.
#' @param col Vector of colors for the different genotypes.
#' @param shift If TRUE, shift the chromosomes so they all start at 0.
#' @param chrwidth Total width of rectangles for each chromosome, as a
#'     fraction of the distance between them.
#' @param ... Additional graphics parameters
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @export
#' @importFrom ggplot2 aes element_rect facet_wrap geom_point geom_rect ggplot scale_y_reverse theme xlab ylab 
#' @importFrom dplyr bind_rows filter mutate rename
#' @importFrom purrr map transpose
#' @importFrom rlang .data
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
#' # inferred genotype at a 28.6 cM on chr 16
#' geno <- maxmarg(probs)
#'
#' ggplot_onegeno(geno, map, shift = TRUE)
#' 
#' ggplot_onegeno(geno, map, ind=1:4)
ggplot_onegeno <-
    function(geno, map, ind=1, chr=NULL, col=NULL,
             shift=FALSE,
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
    
    ggplot_onegeno_internal(geno, map, ind, 
                          col=col, chrwidth=chrwidth, ...)
}

# @param na_col Color for missing segments.
# @param border Color of outer border around chromosome rectangles.
# @param bgcolor Color for background rectangle
ggplot_onegeno_internal <-
  function(geno, map, ind, col=NULL, na_col="white",
           border="black", bgcolor="gray90",
           chrwidth=0.5,
           xlab="Chromosome", ylab="Position (Mbp)",
           ylim=NULL, main="",
           vlines.col="gray80",
           legend.position = "none",
           colors = qtl2::CCcolors,
           ...)
  {
    dots <- list(...)
    if(is.null(ylim)) ylim <- rev(range(unlist(map), na.rm=TRUE))
    
    nchr <- length(map)
    chrwidth <- chrwidth / 2
    
    intervals <- get_geno_intervals(geno, map, ind, chrwidth, colors)
    
    ## initial plot setup
    p <- ggplot2::ggplot(intervals$map) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data$chr, y = .data$lo),
        col = "transparent") +
      ggplot2::ylab(xlab) +
      ggplot2::xlab(ylab) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = unclass(.data$chr) - chrwidth,
          xmax = unclass(.data$chr) + chrwidth, 
          ymin = .data$lo,
          ymax = .data$hi),
        fill = na_col, col = border) +
      ggplot2::theme(
        legend.position = legend.position,
        panel.background = ggplot2::element_rect(fill = bgcolor)) +
      ggplot2::scale_y_reverse() +
      ggplot2::facet_wrap(~ ind)
    
    # grid lines
    p <- ggplot_grid_lines(p, vlines.col=vlines.col, ...)
    
    # color
    max_geno <- max(unlist(geno), na.rm=TRUE)
    if(is.null(col)) {
      if(max_geno <= length(colors)) {
        col <- colors
      }
      else {
        warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
        col <- rep(colors, max_geno)
      }
    }
    else if(max_geno > length(col)) {
      warning("not enough colors; recycling them")
      col <- rep(col, max_geno)
    }
    p <- p + 
      ggplot2::scale_fill_manual(values = col)
    
    # colors by genotype interval
    if(is.null(intervals$right)) {
      p <- p + ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$chrleft,
          xmax = .data$chrright, 
          ymin = .data$lo,
          ymax = .data$hi,
          fill = .data$fill),
        intervals$left) +
        # redraw border
        ggplot2::geom_rect(
          ggplot2::aes(
            xmin = .data$chrleft,
            xmax = .data$chrright, 
            ymin = .data$lo,
            ymax = .data$hi),
          fill = "transparent", col = border)
    } else {
      p <- p + ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$chrleft,
          xmax = .data$chrcode, 
          ymin = .data$lo,
          ymax = .data$hi,
          fill = .data$fill),
        intervals$left) +
        # redraw border
        ggplot2::geom_rect(
          ggplot2::aes(
            xmin = .data$chrleft,
            xmax = .data$chrcode, 
            ymin = .data$lo,
            ymax = .data$hi),
          fill = "transparent", col = border)
      p <- p + ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$chrcode,
          xmax = .data$chrright, 
          ymin = .data$lo,
          ymax = .data$hi,
          fill = .data$fill),
        intervals$right) +
        # redraw border
        ggplot2::geom_rect(
          ggplot2::aes(
            xmin = .data$chrcode,
            xmax = .data$chrright, 
            ymin = .data$lo,
            ymax = .data$hi),
          fill = "transparent", col = border)
    }
    
    # add box just in case
    p + ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = border,
                                           fill=NA))
  }

get_geno_intervals <- function(geno, map, ind = 1, chrwidth = 0.25, colors) {
  # set up genotype intervals
  # for now, geno is reduced in parent function to one interval
  # want to be able to do this for multiple individuals.
  if(is.matrix(geno[[1]])){
    tmpfn <- function(x, ind) 
      geno2intervals(x$geno[ind,, drop = FALSE], x$map)
    intervals <- get_geno_intervals_one(geno, map, ind, tmpfn, chrwidth, colors)
    intervals2 <- NULL
  } else {
    if(!is.array(geno[[1]]) || length(dim(geno[[1]])) != 3 ||
       dim(geno[[1]])[3] != 2)
      stop("geno should be an array individuals x positions x 2 haplotypes")

    tmpfn <- function(x, ind) 
      geno2intervals(x$geno[,,1][ind,, drop = FALSE], x$map)
    intervals <- get_geno_intervals_one(geno, map, ind, tmpfn, chrwidth, colors)
    tmpfn <- function(x, ind) 
      geno2intervals(x$geno[,,2][ind,, drop = FALSE], x$map)
    intervals2 <- get_geno_intervals_one(geno, map, ind, tmpfn, chrwidth, colors)
  }

  map <- 
    dplyr::rename(
      dplyr::mutate(
        dplyr::filter(
          dplyr::bind_rows(
            purrr::map(map, function(x) data.frame(t(range(x)))),
            .id = "chr"),
          .data$chr %in% intervals$chr),
        chr = factor(.data$chr, .data$chr),
        chrcode = as.numeric(unclass(.data$chr)),
        chrleft = .data$chrcode - chrwidth,
        chrright = .data$chrcode + chrwidth),
      lo = .data$X1, hi = .data$X2)
  
  list(map = map, left = intervals, right = intervals2)
}
get_geno_intervals_one <- function(geno, map, ind, tmpfn, chrwidth, colors) {
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
    chr = factor(.data$chr, names(map)[names(map) %in% .data$chr]),
    chrcode = as.numeric(unclass(.data$chr)),
    chrleft = .data$chrcode - chrwidth,
    chrright = .data$chrcode + chrwidth,
    fill = factor(names(colors)[geno], names(colors)))
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

    data.frame(lo = map[c(1, xo_int + 1)],
               hi = map[c(xo_int, length(map))],
               geno = geno[c(xo_int, length(map))])

}
