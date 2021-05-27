#' GGPlot internal routine for ggplot_genes
#'
#' Plot genes at positions
#'
#' @param start,end,strand,rect_top,rect_bottom,colors,space,y,dir_symbol,name,xlim usual parameters
#' @param legend.position,vlines,xlab,ylab,bgcolor,xat hidden parameters
#' @param ... Additional graphics parameters.
#' 
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @importFrom ggplot2 aes coord_cartesian element_blank element_rect geom_text ggplot scale_x_continuous theme xlab ylab
#' geom_rect geom_text 
#' scale_x_continuous 
#' theme element_rect element_blank xlim
ggplot_genes_internal <-
  function(start, end, strand, rect_top, rect_bottom, 
           colors, space, y, dir_symbol, name, xlim,
           ## remaining options are ...
           xlab="Position (Mbp)", ylab="",
           bgcolor="gray92", xat=NULL,
           legend.position = "none",
           vlines=NULL, ...)
  {
    dat <- data.frame(start=start, end=end, strand=strand,
                     rect_top=rect_top, rect_bottom=rect_bottom,
                     colors=colors, name=name, y=y)
    p <- ggplot2::ggplot(dat, 
                         ggplot2::aes(x=end+space, y=y, 
                                      xmin=start,
                                      xmax=end,
                                      ymin=rect_bottom,
                                      ymax=rect_top, 
                                      col=colors,
                                      fill=colors)) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::theme(legend.position = legend.position)
    
    # gray background
    if(!is.null(bgcolor))
      p <- p + ggplot2::theme(panel.background = 
                                ggplot2::element_rect(fill = bgcolor,
                                                      color = "black"))
    # axis
    if(is.null(xat)) xat <- pretty(xlim)
    if(length(xat) > 1 || !is.na(xat)) {
      xlim <- range(xat)
      p <- p + 
        ggplot2::scale_x_continuous(breaks = xat) +
        ggplot2::coord_cartesian(xlim = xlim)
    }  

    # vertical lines
    if(length(vlines)==1 && is.na(vlines)) {
      # if vlines==NA, skip lines
      p <- p + 
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank())
    }

    p <- p + 
      ggplot2::geom_rect() +
      # gene symbol
      ggplot2::geom_text(mapping = 
                           ggplot2::aes(x = end + space,
                                        y = y,
                                        label = paste0("'", name, "'", dir_symbol),
                                        hjust = 0,
                                        vjust = 0.5),
                         parse = TRUE)

    # add black box
    p +
      ggplot2::theme(panel.border = 
                       ggplot2::element_rect(colour = "black",
                                             fill=NA))
  }
