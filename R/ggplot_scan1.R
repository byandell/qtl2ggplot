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
#' @param box Include box.
#'
#' @param ... Additional graphics parameters.
#' 
ggplot_scan1 <-
  function(map, lod, gap, group = NULL,
           bgcolor, altbgcolor,
           lwd=1, pch=1, cex=0.5, col=NULL, xlab=NULL, ylab="LOD score",
           xaxt = "y", yaxt = "y",
           palette = "Dark2",
           xlim=NULL, ylim=NULL, main="",
           hlines=NULL, vlines=NULL,
           legend.position = 
             ifelse(ncol(lod) == 1, "none", "right"),
           legend.title="pheno",
           lines=TRUE, points=!lines,
           box=TRUE,
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
    df <- data.frame(xpos=xpos, chr=chr, lod) %>%
      gather(pheno, lod, -xpos, -chr) %>%
      mutate(pheno = as.character(pheno))

    # Use group if provided and only one pheno
    if(!is.null(group)) {
      if(length(group) == nrow(df) & ncol(lod) == 1)
        df$pheno <- factor(group)
    } else {
      df$pheno <- factor(df$pheno)
      levels(df$pheno) <- dimnames(lod)[[2]]
    }
    df <- df %>%
      mutate(group = paste(chr, pheno, sep = "_"))

    # make ggplot aesthetic with limits and labels
    p <- ggplot(df, aes(x=xpos,y=lod,col=pheno,group=group)) +
      ylim(ylim) +
      xlab(xlab) +
      ylab(ylab)

    ## Add lines and/or points.
    if(lines) {
      p <- p + geom_line(size = lwd)
    }
    if(points) {
      p <- p + geom_point(shape = pch,
                          size = cex)
    }

    if(is.null(col)) {
      if(is.null(palette)) palette <- "Dark2"
      p <- p + scale_color_brewer(name = legend.title,
                                  palette = palette)
    } else {
      col <- rep(col, length = length(unique(df$pheno)))
      names(col) <- NULL
      p <- p + scale_color_manual(name = legend.title,
                                  values = col)
    }


    # add legend if requested
    p <- p +
      theme(legend.position = legend.position)

    # add background rectangles
    if(!is.null(bgcolor))
      p <- p + theme(panel.background = element_rect(fill = bgcolor))
    if(!is.null(altbgcolor) && !onechr) {
      df_rect <- data.frame(xmin=chrbound[1,], xmax=chrbound[2,],
                            ymin=ylim[1], ymax=ylim[2])
      df_rect <- df_rect[seq(2, ncol(chrbound), by=2),]
      # Not sure why color,x,y needed in geom_rect
      p <- p + geom_rect(mapping = aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                                       # need to unmap following
                                       color=NULL,x=NULL,y=NULL,group=NULL),
                         data = df_rect,
                         fill = altbgcolor, col = altbgcolor)
    }

    # include axis labels?
    if(yaxt == "n") {
      p <- p + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
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
      p <- p + scale_x_continuous(breaks = loc,
                                  labels = names(map),
                                  lim = xlim)
    }

    # remove y axis?
    if(yaxt == "n") {
      p <- p + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    }
    # grid lines
    if((length(vlines)==1 && is.na(vlines)) | !onechr) {
      # if vlines==NA (or mult chr), skip lines
      p <- p + theme(panel.grid.major.x=element_blank(),
                     panel.grid.minor.x=element_blank())
    }
    if((length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
      p <- p + theme(panel.grid.major.y=element_blank(),
                     panel.grid.minor.y=element_blank())
    }
    # add box just in case
    if(box) {
      p <- p +
        theme(panel.border = element_rect(colour = "black",
                                          fill=NA))
    }
    # add main as title if provided
    # or use name from lod if only one column
    if(!is.logical(main)) {
      title <- main
      main <- TRUE
    }
    if(main) {
      if(title == "") {
        # create title from pheno name if only 1
        pheno_names <- levels(df$pheno)
        if(length(pheno_names) == 1) {
          p <- p +
            ggtitle(pheno_names) +
            theme(legend.position = "none")
        }
      } else {
        p <- p + ggtitle(title)
      }
    }

    p
  }
