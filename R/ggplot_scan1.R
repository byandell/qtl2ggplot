#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param map Map of pseudomarker locations.
#' @param lod Matrix of lod (or other) values.
#' @param gap Gap between chromosomes.
#' @param col colors for points or lines, with labels.
#' @param pattern Use to group values for plotting (default = \code{NULL}); typically provided by \code{\link{plot_snpasso}} internal routine.
#' @param facet Plot facets if multiple phenotypes and group provided (default = \code{NULL}).
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}.
#'
#' @param ... Additional graphics parameters.
#' 
#' @param bgcolor Background color for the plot.
#' @param altbgcolor Background color for alternate chromosomes.
#' @param lwd,pch,cex,col,xlim,ylim,xaxt,yaxt base plot parameters (coverted to ggplot use)
#' @param palette color palette for points and lines
#' @param xlab,ylab,main Titles for x,y and plot.
#' @param hlines,vlines Horizontal and vertical lines.
#' @param legend.position,legend.title Legend theme setting.
#' @param lines,points Include lines and/or points.
#' 
#' @importFrom ggplot2 ggplot aes xlim ylim xlab ylab ggtitle 
#' geom_line geom_point theme geom_rect facet_wrap
#' scale_x_continuous
#' theme element_rect element_blank
#' @importFrom tidyr gather
#' @importFrom dplyr mutate rename
ggplot_scan1 <-
  function(map, lod, gap,
           col=NULL, 
           shape=NULL,
           pattern = NULL, facet = NULL,
           patterns = c("none","all","hilit"),
           ...)
  {
    patterns <- match.arg(patterns)
    scan1ggdata <- make_scan1ggdata(map, lod, gap, col, pattern, shape,
                                    facet, patterns)
    
    ggplot_scan1_internal(map, gap, col, shape, scan1ggdata, facet, ...)
  }

make_scan1ggdata <- function(map, lod, gap, col, pattern, shape,
                             facet, patterns) {
  # set up chr and xpos with gap.
  xpos <- map_to_xpos(map, gap)
  chr <- rep(names(map), sapply(map, length))
  
  # make data frame for ggplot
  scan1ggdata <- data.frame(xpos=xpos, chr=chr, lod,
                            check.names = FALSE)
  scan1ggdata <- tidyr::gather(scan1ggdata, pheno, lod, -xpos, -chr)
  # make sure order of pheno is preserved.
  scan1ggdata <- dplyr::mutate(scan1ggdata, 
                               pheno = ordered(pheno, levels = unique(pheno)))
  
  ## facet if more than one pheno or set by user.
  if(ncol(lod) > 1 & !is.null(pattern)) {
    # If facet is not NULL, pattern  column of scan1ggdata is used to facet.
    # That column is either pheno or pattern, set in color_patterns_pheno.
    if(is.null(facet))
      facet <- "pheno"
  }
  ## Set up col, group and (optional) facet in scan1ggdata.
  ## Column pheno becomes either col or facet
  color_patterns_pheno(scan1ggdata,
                       lod, pattern, col, shape,
                       patterns, facet)
}

ggplot_scan1_internal <-
  function(map, gap, col, shape, scan1ggdata, facet,
           bgcolor, altbgcolor,
           lwd=1, 
           pch = c(SNP=96,indel=23,INS=25,DEL=24,INV=22), 
           cex=1, 
           xlab=NULL, ylab="LOD score",
           xaxt = "y", yaxt = "y",
           palette = "Dark2",
           xlim=NULL, ylim=NULL, main=FALSE,
           hlines=NULL, vlines=NULL,
           legend.position = 
             ifelse(length(levels(scan1ggdata$color)) == 1, "none", "right"),
           legend.title="pheno",
           lines=TRUE, points=!lines,
           ...)
  {
    
    # Extra arguments
    onechr <- (length(map)==1) # single chromosome

    chrbound <- map_to_boundaries(map, gap)

    if(is.null(ylim))
      ylim <- c(0, max(scan1ggdata$lod, na.rm=TRUE)*1.02)

    if(is.null(xlim)) {
      xlim <- range(scan1ggdata$xpos, na.rm=TRUE)
      if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
    }

    if(is.null(xlab)) {
      if(onechr) {
        if(names(map) == " ") xlab <- "Position"
        else xlab <- paste("Chr", names(map), "position")
      }
      else xlab <- "Chromosome"
    }

    ## filter data so only using what we will plot.
    scan1ggdata <- dplyr::filter(scan1ggdata, 
                                 lod >= ylim[1] & lod <= ylim[2])
    
    # make ggplot aesthetic with limits and labels
    p <- ggplot2::ggplot(scan1ggdata, 
                         ggplot2::aes(x = xpos, y = lod,
                                      col = color, 
                                      shape = shape,
                                      group = group)) +
      ggplot2::ylim(ylim) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
    
    # Facets (if multiple phenotypes and groups).
    if(!is.null(facet)) {
      p <- p + ggplot2::facet_wrap(~ facets)
    }

    # color palette, point shapes and legend titles
    col_shape <- color_patterns_get(scan1ggdata, col, palette, shape)
    p <- p +
      ggplot2::scale_color_manual(name = legend.title,
                                  values = col_shape$colors)
    if(length(col_shape$shapes) > 1)
    p <- p +
      ggplot2::scale_shape_manual(name = "SV Type",
                                  labels = names(col_shape$shapes),
                                  values = col_shape$shapes)
    
    # add legend if requested
    p <- p +
      ggplot2::theme(legend.position = legend.position)

    # add background rectangles
    if(!is.null(bgcolor)) {
      p <- p + 
        ggplot2::theme(panel.background = 
                         ggplot2::element_rect(fill = bgcolor))
    }
    if(!is.null(altbgcolor) && !onechr) {
      df_rect <- data.frame(xmin=chrbound[1,], xmax=chrbound[2,],
                            ymin=ylim[1], ymax=ylim[2])
      df_rect <- df_rect[seq(2, ncol(chrbound), by=2),]
      p <- p + 
        ggplot2::geom_rect(mapping = aes(xmin=xmin, 
                                         xmax=xmax, 
                                         ymin=ymin, 
                                         ymax=ymax),
                           inherit.aes = FALSE,
                           data = df_rect,
                           fill = altbgcolor, 
                           col = altbgcolor)
    }

    # include axis labels?
    if(yaxt == "n") {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank())
    }
    if(onechr) {
      if(xaxt == "n") {
        p <- p + theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
      }
      p <- p + xlim(xlim)
    } else {
      # x axis for multiple chromosomes
      loc <- colMeans(chrbound)
      p <- p + 
        ggplot2::scale_x_continuous(breaks = loc,
                                    labels = names(map),
                                    lim = xlim)
    }

    # remove y axis?
    if(yaxt == "n") {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank())
    }
    # grid lines
    if((length(vlines)==1 && is.na(vlines)) | !onechr) {
      # if vlines==NA (or mult chr), skip lines
      p <- p + 
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank())
    }
    if((length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
      p <- p + 
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank())
    }
    # add box just in case
    p <- p +
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black",
                                        fill=NA))

    # add main as title if provided
    # or use name from lod if only one column
    if(!is.logical(main)) {
      title <- main
      main <- TRUE
    }
    if(main) {
      if(title == "") {
        # create title from pheno name if only 1
        group_names <- levels(df$pheno)
        if(length(group_names) == 1) {
          p <- p +
            ggplot2::ggtitle(group_names) +
            ggplot2::theme(legend.position = "none")
        }
      } else {
        p <- p + 
          ggplot2::ggtitle(title)
      }
    }
    
    ## Add lines and/or points.
    if(lines) {
      p <- p + ggplot2::geom_line(size = lwd)
    }
    if(points) {
      p <- p + ggplot2::geom_point(shape = pch,
                                   size = cex, 
                                   fill = "grey40")
    }
    
    p
  }
