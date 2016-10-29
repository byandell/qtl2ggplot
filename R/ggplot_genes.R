ggplot_genes <-
  function(start, end, strand, rect_top, rect_bottom, 
           colors, space, y, name, 
           ## remaining options are ...
           xlab="Position (Mbp)", xaxs="i",
           bgcolor="gray92", xat=NULL,
           mgp=c(0,0.2,0),
           vlines=NULL, vlines.col="white",
           vlines.lwd=1, vlines.lty=1)
  {
    ## Look at doqtl2::plot.feature_tbl
    
    df <- data.frame(start=start, end=end, strand=strand,
                     rect_top=rect_top, rect_bottom=rect_bottom,
                     colors=colors, name=name)
    p <- ggplot(df, aes(xmin=start, xmax=end,
                        ymin=rect_bottom, ymax=rect_top, 
                        col=colors)) +
      xlab(xlab) +
      ylab(ylab)

    # gray background
    u <- par("usr")
    rect(u[1], u[3], u[2], u[4], col=bgcolor, border="black")
    
    # axis
    if(is.null(xat)) xat <- pretty(xlim)
    if(length(xat) > 1 || !is.na(xat))
      axis(side=1, at=xat, mgp=mgp, tick=FALSE)
    
    # vertical lines
    if(is.null(vlines)) vlines <- xat
    if(length(vlines) > 1 || !is.na(vlines))
      abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
    
    for(i in seq(along=start)) {
      rect(start[i], rect_top[i],
           end[i],   rect_bottom[i],
           col=colors[i], border=colors[i],
           lend=1, ljoin=1)
      text(end[i] + space, y[i],
           name[i], adj=c(0, 0.5), col=colors[i],
           cex=text_cex)
      if(!is.na(strand[i]) && (strand[i] == "+" || strand[i] == '-'))
        text(end[i] + space + strwidth(name[i], cex=text_cex), y[i],
             dir_symbol[i], adj=c(0, 0.5), col=colors[i],
             cex=text_cex)
    }

    box()
  }
