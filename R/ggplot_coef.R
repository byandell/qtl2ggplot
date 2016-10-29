ggplot_coef <-
  function(x, columns, ylim, col, add, gap, bgcolor, altbgcolor, ylab, 
           legend.position = "right", legend.title = "geno", ...)
  {
    plot_scan1(x, lodcolumn=columns, ylim=ylim, col=col, add=add,
               gap=gap, bgcolor=bgcolor, altbgcolor=altbgcolor,
               ylab=ylab, 
               legend.position=legend.position, 
               legend.title = legend.title,
               ...)
  }
