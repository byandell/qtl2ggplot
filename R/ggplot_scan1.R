ggplot_scan1 <-
  function(map, lod, add=FALSE, gap,
           bgcolor, altbgcolor,
           lwd=2, col=NULL, xlab=NULL, ylab="LOD score",
           xlim=NULL, ylim=NULL, main="",
           hlines=NULL, vlines=NULL,
           legend.position="none",
           legend.title="pheno",
           lines=TRUE, points=FALSE,
           ...)
  {
    # Extra arguments
    dots <- list(...)
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
    chr <- rep(names(map), sapply(map, length))
    df <- data.frame(xpos=xpos, chr=chr, lod) %>%
      gather(pheno,lod,-xpos,-chr)
    # Use group if provided and only one pheno 
    if(!is.null(dots$group)) {
      if(length(dots$group) == nrow(df) & ncol(lod) == 1)
        df$pheno <- dots$group
    }
    df <- df %>%
      mutate(group = paste(df$chr, df$pheno, sep = "_"))
    
    # make ggplot aesthetic with limits and labels
    p <- ggplot(df, aes(x=xpos,y=lod,col=pheno,group=group)) +
      xlim(xlim) +
      ylim(ylim) +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(main)
    
    if(is.null(col)) {
      if(is.null(dots$palette)) dots$palette <- "Dark2"
      p <- p + scale_color_brewer(name = legend.title,
                                  palette = dots$palette)
    } else {
#      col <- rep(col, length = ncol(lod))
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
                                       color=altbgcolor,x=xmin,y=ymin,group="1"),
                         data = df_rect,
                         fill = altbgcolor, col = altbgcolor)
    }
    
    # include axis labels?
    if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
    if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")
    if(dots$yaxt == "n") {
      p <- p + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    }
    if(onechr) {
      if(dots$xaxt == "n") {
        p <- p + theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
      }
    } else {
      # x axis for multiple chromosomes
      loc <- colMeans(chrbound)
      p <- p + scale_x_continuous(breaks = loc,
                                  labels = names(map))
    }
    
    # add y axis unless par(yaxt="n")
    if(dots$yaxt == "n") {
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
    p <- p + geom_rect(aes(xmin=xlim[1], xmax=xlim[2], 
                           ymin=ylim[1], ymax=ylim[2]),
                           fill = NA, col = "black")
    if(lines)
      p <- p + geom_line()
    if(points) {
      if(is.null(dots$pch)) dots$pch <- 1
      if(is.null(dots$cex)) dots$cex <- 0.5
      p <- p + geom_point(shape = dots$pch,
                          size = dots$cex)
    }
    p
  }
