ggplot_genes <-
  function(start, end, strand, rect_top, rect_bottom, 
           colors, space, y, dir_symbol, name, xlim,
           ## remaining options are ...
           xlab="Position (Mbp)", ylab="",
           bgcolor="gray92", xat=NULL,
           legend.position = "none",
           vlines=NULL, vlines.col="white",
           vlines.lwd=1, vlines.lty=1)
  {
    ## Look at doqtl2::plot.feature_tbl
    
    df <- data.frame(start=start, end=end, strand=strand,
                     rect_top=rect_top, rect_bottom=rect_bottom,
                     colors=colors, name=name, y=y)
    p <- ggplot(df, aes(x=end+space, y=y, xmin=start, xmax=end,
                        ymin=rect_bottom, ymax=rect_top, 
                        col=colors, fill=colors)) +
      xlab(xlab) +
      ylab(ylab) +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      theme(legend.position = legend.position)

    # gray background
    if(!is.null(bgcolor))
      p <- p + theme(panel.background = element_rect(fill = bgcolor,
                                                     color = "black"))
    # axis
    if(is.null(xat)) xat <- pretty(xlim)
    if(length(xat) > 1 || !is.na(xat))
      p <- p + scale_x_continuous(breaks = xat)

    # vertical lines
    if(length(vlines)==1 && is.na(vlines)) {
      # if vlines==NA, skip lines
      p <- p + theme(panel.grid.major.x=element_blank(),
                     panel.grid.minor.x=element_blank())
    }

    p <- p + 
      geom_rect() +
      # gene symbol
      geom_text(mapping = aes(x = end + space,
                              y = y,
                              label = paste0("'", name, "'", dir_symbol),
                              hjust = 0,
                              vjust = 0.5),
                parse = TRUE)

    # add black box
    p + geom_rect(aes(xmin=xlim[1], xmax=xlim[2], 
                      ymin=0, ymax=1),
                      fill = NA, col = "black")
  }
