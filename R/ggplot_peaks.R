#' Plot QTL peak locations
#'
#' Plot QTL peak locations (possibly with intervals) for multiple traits.
#'
#' @param peaks Data frame such as that produced by
#'     \code{\link[qtl2]{find_peaks}}) containing columns
#'     \code{chr}, \code{pos}, \code{lodindex}, and \code{lodcolumn}.
#'     May also contain columns \code{ci_lo} and \code{ci_hi}, in
#'     which case intervals will be plotted.
#' @param map Marker map, used to get chromosome lengths (and start
#'     and end positions).
#' @param chr Selected chromosomes to plot; a vector of character
#'     strings.
#' @param tick_height Height of tick marks at the peaks (a number between 0 and 1).
#' @param bgcolor Background color for the plot.
#' @param altbgcolor Background color for alternate chromosomes.
#' @param gap Gap between chromosomes.
#' @param ... Additional graphics parameters
#'
#' @seealso \code{\link[qtl2]{find_peaks}}
#'
#' @export
#' @importFrom ggplot2 aes element_blank facet_grid geom_segment 
#' ggplot ggtitle theme xlab ylab
#' @importFrom grid unit
#' @importFrom rlang .data
#' 
#' @return None.
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
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # find peaks above lod=3.5 (and calculate 1.5-LOD support intervals)
#' peaks <- find_peaks(out, map, threshold=3.5, drop=1.5)
#'
#' # color peaks above 6 red; only show chromosomes with peaks
#' plot_peaks(peaks, map)
#' peaks$col <- (peaks$lod > 6)
#' 
#' ggplot_peaks(peaks, map[names(map) %in% peaks$chr], col = c("blue","red"),
#'            legend.title = "LOD > 6")

ggplot_peaks <- 
  function(peaks, map, chr=NULL, tick_height = 0.3,
           gap=25, bgcolor="gray90", altbgcolor="gray85",
           ...) 
  {
    if(is.null(map)) stop("map is NULL")
    
    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names
    
    if(!is.null(chr)) {
      chri <- match(chr, names(map))
      if(any(is.na(chri)))
        stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
      map <- map[chri]
      peaks <- peaks[peaks$chr %in% names(map),,drop=FALSE]
    }
    
    if(nrow(peaks) == 0)
      stop("There are no peaks on the selected chromosomes")
    
    ggplot_peaks_internal(peaks=peaks, map=map, tick_height=tick_height,
                        gap=gap, bgcolor=bgcolor, algbgcolor=altbgcolor, ...)
    
  }

ggplot_peaks_internal <-
  function(peaks, map, tick_height,
           gap, bgcolor, algbgcolor,
           lwd=2, col="slateblue", xlab=NULL, ylab="",
           xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
           main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
           mgp=NULL, las=1, lend=1, ljoin=1,
           hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
           vlines=NULL, vlines.col="white", vlines.lwd=1, vlines.lty=1,
           point_size=0, vbars=TRUE,
           xaxt=ifelse(onechr, "y", "n"),
           yaxt="y",
           legend.title = "",
           legend.position = "right",
           ...)
  {
    dots <- list(...)
    onechr <- (length(map)==1) # single chromosome
    
    # make chr into factor
    chrs <- names(map)
    peaks$chr <- factor(peaks$chr, chrs)
    peaks$lodcolumn <- factor(peaks$lodcolumn)
    
    # color
    if(is.null(peaks$col)) {
      peaks$col <- "all"
      legend.position <- "none"
    }
    
    # get map lengths by chr.
    mapl <- t(sapply(map[chrs], range))
    mapl <- data.frame(mapl)
    colnames(mapl) <- c("lo","hi")
    mapl$chr <- factor(chrs, chrs)
    
    if(is.null(xlab)) {
      if(onechr) {
        if(names(map) == " ") xlab <- "Position"
        else xlab <- paste("Chr", names(map), "position")
      }
      else xlab <- "Chromosome"
    }
    
    p <- ggplot2::ggplot(peaks) +
      ggplot2::aes(x = .data$pos, y = .data$lodcolumn, group = .data$chr, col = col) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
    
    # Add color and legend.
    p <- p +
      ggplot2::theme(legend.position = legend.position) +
      ggplot2::scale_color_manual(name = legend.title,
                                  values = col)
    
    if(main != "") {
      p <- p +
        ggplot2::ggtitle(main)
    }
    
    if(!onechr) {
      p <- p +
        ggplot2::facet_grid(~ .data$chr, scales = "free_x", space = "free_x") +
        ggplot2::theme(panel.spacing = grid::unit(gap / 10000, "npc"))
    }
    
    p <- p +   
      ggplot2::geom_point(size = point_size)
    
    # set up horizontal axis to match data.
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(x = .data$lo, xend = .data$hi,
                                         y = 1, yend = 1),
                            data = mapl, col = "transparent",
                            inherit.aes = FALSE)
    
    # include axis labels?
    if(yaxt == "n") {
      p <- p + 
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank())
    }
    # X axis
    if(xaxt == "n") {
      p <- p + 
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank())
    }
    
    # grid lines
    if((length(vlines)==1 && is.na(vlines)) | !onechr) { # if vlines==NA (or mult chr), skip lines
      p <- p +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank())
    }
    if((length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
      p <- p +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.y = ggplot2::element_blank())
    }
    
    if(vbars) {
      tick_height <- tick_height / 2
      p <- p +
        ggplot2::geom_segment(
          ggplot2::aes(xend = .data$pos,
                       y    = unclass(.data$lodcolumn) - tick_height,
                       yend = unclass(.data$lodcolumn) + tick_height))
    }
    if("ci_lo" %in% names(peaks) && "ci_hi" %in% names(peaks)) {
      p <- p +
        ggplot2::geom_segment(
          ggplot2::aes(x = .data$ci_lo, xend = .data$ci_hi,
                       y = .data$lodcolumn, yend = .data$lodcolumn))
    }
    # Add box for each chr.
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black",
                                             fill=NA))
    
    p
  }
