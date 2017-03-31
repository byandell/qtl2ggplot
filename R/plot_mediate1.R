#' Plot mediate1 object
#' 
#' Plot mediation by expression phenotypes in region.
#' 
#' @param x Object of class \code{\link[qtl2pattern]{allele1}}.
#' @param ... Other parameters ignored.
#' 
#' @export
#' @importFrom ggplot2 aes autoplot geom_hline geom_point geom_vline ggplot ggtitle
#' 
plot_mediate1 <- function(x, ...) {
  pos <- attr(x, "pos")
  lod <- attr(x, "lod")
  ggplot2::ggplot(x, ggplot2::aes(x=start, y=lod_t_m, symbol = symbol)) +
    ggplot2::geom_hline(yintercept = lod, col = "darkgray") +
    ggplot2::geom_vline(xintercept = pos, col = "darkgray") +
    ggplot2::geom_point(shape = 1) +
    ggplot2::ggtitle(attr(x, "pheno"))
}
#' @export
plot.mediate1 <- function(x, ...)
  ggplot2::autoplot(x, ...)
#' @export
autoplot.mediate1 <- function(x, ...)
  plot_mediate1(x, ...)
