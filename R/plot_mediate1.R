#' Plot mediate1 object
#' 
#' Plot mediation by expression phenotypes in region.
#' 
#' @param x Object of class \code{\link[qtl2pattern]{allele1}}.
#' @param frame If \code{TRUE}, enable frames for \code{\link[plotly]{ggplotly}}
#' @param ... Other parameters ignored.
#' 
#' @export
#' @importFrom ggplot2 aes autoplot geom_hline geom_point geom_vline ggplot ggtitle
#' 
plot_mediate1 <- function(x, frame = FALSE, ...) {
  
  # Frame does not work properly with plotly::ggplotly.
  # Explore putting frame in geom_point or ...?
  
  pos_raw <- attr(x, "pos")
  lod_raw <- attr(x, "lod")

  scans <- as.data.frame(x$scan)
  scans$pos <- x$map
  scans <- tidyr::gather(scans, id, lod_t_m, -pos)
  if(!frame)
    scans <- dplyr::filter(scans,
                           pos == pos_raw)
  scans <- dplyr::full_join(x$annot, scans, by = "id")

  p <- 
    ggplot2::ggplot(scans, 
      ggplot2::aes(x=start, y=lod_t_m, symbol = symbol)) +
    ggplot2::geom_vline(xintercept = pos_raw, col = "darkgray") +
    ggplot2::geom_hline(yintercept = lod_raw, col = "darkgray") +
    ggplot2::ggtitle(attr(x, "pheno"))
  
  if(frame) {
    p <- p +
      ggplot2::geom_point(shape = 1, size = 3, 
                          ggplot2::aes(frame = pos))
  } else {
    p <- p +
      ggplot2::geom_point(shape = 1, size = 3)
  }
  p
}
#' @export
plot.mediate1 <- function(x, ...)
  ggplot2::autoplot(x, ...)
#' @export
autoplot.mediate1 <- function(x, ...)
  plot_mediate1(x, ...)
