#' @export
#' @importFrom ggplot2 aes element_blank 
#' facet_grid geom_text ggplot scale_x_continuous theme
#' 
plot.allele1 <- function(x, ...) {
  autoplot.allele1(x, ...)
}
autoplot.allele1 <- function(x, ...) {
  ggplot2::ggplot(x, 
                  ggplot2::aes(x=1, y=effect, label=allele, col = probe)) +
    ggplot2::geom_text(size=3, position = ggplot2::position_jitter()) +
    ggplot2::facet_grid(~source, scales = "free") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand=c(0,0.1))
} 
