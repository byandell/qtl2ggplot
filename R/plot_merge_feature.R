#' Plot of merge_feature object
#'
#' @param x of class \code{merge_feature}
#' @param pheno name of phenotype to be plotted
#' @param plot_by element to plot by (one of \code{c("pattern","consequence")})
#' @param ... other arguments not used
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 aes facet_wrap geom_jitter ggplot ggtitle xlab ylab
#'
plot_merge_feature <- function(x, pheno, plot_by=c("pattern","consequence"), ...) {
  plot_by <- match.arg(plot_by)
  x$lod <- x[[pheno]]
  x <- dplyr::filter(
    dplyr::mutate(x, pattern = sdp_to_pattern(sdp)),
    !is.na(lod))
  switch(plot_by,
         pattern = {
           ggplot2::ggplot(x,
                           ggplot2::aes(x=pos_Mbp,y=lod,col=pattern)) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~snp_type, scale = "free") +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Consequence for", pheno))
         },
         consequence = {
           ggplot2::ggplot(x,
                           ggplot2::aes(x=pos_Mbp,y=lod,col=snp_type)) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~pattern, scale = "free") +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Allele Pattern for", pheno))
         })
}

#' @method autoplot merge_feature
#' @export
#' @export autoplot.merge_feature
#' @rdname plot_merge_feature
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.merge_feature <- function(x, ...)
  plot_merge_feature(x, ...)

#' @method plot merge_feature
#' @export
#' @export plot.merge_feature
#' @rdname plot_merge_feature
#' 
plot.merge_feature <- function(x, ...)
  autoplot.merge_feature(x, ...)
