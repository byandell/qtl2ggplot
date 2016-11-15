#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param map Map of pseudomarker locations.
#' @param lod Matrix of lod (or other) values.
#' @param gap Gap between chromosomes.
#' @param bgcolor Background color for the plot.
#' @param altbgcolor Background color for alternate chromosomes.
#' @param lwd,pch,cex,col,xlim,ylim,xaxt,yaxt base plot parameters (coverted to ggplot use)
#' @param palette color palette for points and lines
#' @param xlab,ylab,main Titles for x,y and plot.
#' @param hlines,vlines Horizontal and vertical lines.
#' @param legend.position,legend.title Legend theme setting.
#' @param lines,points Include lines and/or points.
#' @param group Use to group values for plotting (default = \code{NULL}).
#' @param facet Plot facets if multiple phenotypes and group provided.
#'
#' @param ... Additional graphics parameters.
#' 
#' @importFrom ggplot2 ggplot aes xlim ylim xlab ylab ggtitle 
#' geom_line geom_point theme geom_rect facet_wrap
#' scale_color_brewer scale_color_manual scale_x_continuous
#' theme element_rect element_blank
#' @importFrom tidyr gather
#' @importFrom dplyr mutate filter rename full_join
ggplot_scan1 <-
  function(map, lod, gap,
           bgcolor, altbgcolor,
           lwd=1, pch=1, cex=0.5, col=NULL, xlab=NULL, ylab="LOD score",
           xaxt = "y", yaxt = "y",
           palette = "Dark2",
           xlim=NULL, ylim=NULL, main=FALSE,
           hlines=NULL, vlines=NULL,
           legend.position = 
             ifelse(ncol(lod) == 1, "none", "right"),
           legend.title="pheno",
           lines=TRUE, points=!lines,
           group = NULL, facet = NULL,
           ...)
  {
    # Extra arguments
    onechr <- (length(map)==1) # single chromosome

    xpos <- map_to_xpos(map, gap)
    chrbound <- map_to_boundaries(map, gap)

    if(is.null(ylim))
      ylim <- c(0, max(lod, na.rm=TRUE)*1.02)

    if(is.null(xlim)) {
      xlim <- range(xpos, na.rm=TRUE)
      if(!onechr) xlim <- xlim + c(-gap/2, gap/2)
    }

    if(is.null(xlab)) {
      if(onechr) {
        if(names(map) == " ") xlab <- "Position"
        else xlab <- paste("Chr", names(map), "position")
      }
      else xlab <- "Chromosome"
    }

    # make data frame for ggplot
    # make sure order of pheno is preserved.
    chr <- rep(names(map), sapply(map, length))
    scan1ggdata <- data.frame(xpos=xpos, chr=chr, lod)
    scan1ggdata <- tidyr::gather(scan1ggdata, pheno, lod, -xpos, -chr)
    scan1ggdata <- dplyr::mutate(scan1ggdata, pheno = as.character(pheno))

    # Use group if provided and only one pheno
    if(!is.null(group)) {
      # If provided, group has to be same size as lod.
      stopifnot(length(group) == length(lod))
      if(ncol(lod) == 1) {
        scan1ggdata$pheno <- factor(group)
      } else {
        # Set up facet as either pheno or group.
        group <- as.data.frame(matrix(group, nrow(lod), ncol(lod)))
        names(group) <- names(lod)
        group <- tidyr::gather(group, pheno, lod)
        if(is.null(facet))
          facet <- "pheno"
        if(facet == "pheno") {
          scan1ggdata <- dplyr::rename(scan1ggdata, group = pheno)
          scan1ggdata$pheno <- group$lod
        } else {
          scan1ggdata$group <- group$lod
        }
      }
    } else {
      scan1ggdata$pheno <- factor(scan1ggdata$pheno)
      levels(scan1ggdata$pheno) <- dimnames(lod)[[2]]
    }
    
    # Set up colors if provided.
    if(!is.null(col)) {
      col <- rep(col, length = length(unique(scan1ggdata$pheno)))
      names(col) <- NULL
    }
    ## If no group, and facet provided, set up group using pheno and maybe col.
    if(!is.null(facet) && is.null(scan1ggdata$group)) {
      if(is.null(col))
        scan1ggdata$group <- scan1ggdata$pheno
      else
        scan1ggdata$group <- col[scan1ggdata$pheno]
    }
    
    ## chr_pheno makes chr and pheno combination distinct for plotting.
    scan1ggdata <- dplyr::mutate(scan1ggdata,
                        chr_pheno = paste(chr, pheno, sep = "_"))
    ## filter data so only using what we will plot
    scan1ggdata <- dplyr::filter(scan1ggdata, lod >= ylim[1] & lod <= ylim[2])

    # make ggplot aesthetic with limits and labels
    p <- ggplot2::ggplot(scan1ggdata, 
                         ggplot2::aes(x = xpos, y = lod,
                                      col = pheno, 
                                      group = chr_pheno)) +
      ggplot2::ylim(ylim) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)

    ## Add lines and/or points.
    if(lines) {
      p <- p + ggplot2::geom_line(size = lwd)
    }
    if(points) {
      p <- p + ggplot2::geom_point(shape = pch,
                          size = cex)
    }
    
    ## Facets
    if(!is.null(facet)) {
      p <- p + ggplot2::facet_wrap(~group)
    }

    if(is.null(col)) {
      if(is.null(palette)) palette <- "Dark2"
      p <- p + 
        ggplot2::scale_color_brewer(name = legend.title,
                                    palette = palette)
    } else {
      p <- p + 
        ggplot2::scale_color_manual(name = legend.title,
                                    values = col)
    }


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
                                         ymax=ymax,
                                         # need to unmap following
                                         color=NULL,
                                         x=NULL,
                                         y=NULL,
                                         group=NULL),
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

    p
  }
