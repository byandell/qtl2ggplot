ggplot_grid_lines <- function(p, onechr = TRUE,
                              hlines=NULL, 
                              vlines=NULL, 
                              hlines.col="white", 
                              hlines.lwd=1, 
                              hlines.lty=1,
                              vlines.col="white", 
                              vlines.lwd=1, 
                              vlines.lty=1,
                              ...) {
  # grid lines
  if((length(vlines)==1 && is.na(vlines)) | !onechr) {
    # if vlines==NA (or mult chr), skip lines
    p <- p +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank())
  } else {
    if(!missing(vlines.col) | !missing(vlines.lwd) | !missing(vlines.lty)) {
      p <- p +
        ggplot2::theme(
          panel.grid.major.x = 
            ggplot2::element_line(
              color = vlines.col,
              size = vlines.lwd / 2,
              linetype = vlines.lty))
    }
  }
  if((length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
    p <- p +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank())
  } else {
    if(!missing(hlines.col) | !missing(hlines.lwd) | !missing(hlines.lty)) {
      p <- p +
        ggplot2::theme(
          panel.grid.major.y = 
            ggplot2::element_line(
              color = hlines.col,
              size = hlines.lwd / 2,
              linetype = hlines.lty))
    }
  }
  p
}