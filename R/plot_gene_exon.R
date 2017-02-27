#' Plot of exons for a gene with SNPs
#'
#' Uses \code{\link{gene_plot}} to plot genes, exons, mRNA with SNPs.
#'
#' @param exon_tbl tbl of feature information from \code{\link{get_mgi_features}}
#' @param top_snps_tbl table from \code{\link[qtl2scan]{top_snps}}
#' @param plot_now plot now if TRUE
#' @param ... arguments passed along to \code{\link{gene_plot}}
#'
#' @return list of ggplots (see \code{\link{gene_plot}})
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @export
#' 
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom ggplot2 ggtitle
#' 
plot_gene_exon <- function(gene_exon, top_snps_tbl=NULL, plot_now=TRUE,
                           genes = unique(gene_exon$gene), ...) {
  p <- list()
  ## Reduce to max lod per SNP
  if(!is.null(top_snps_tbl)) {
    top_snps_tbl <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(top_snps_tbl, snp_id,pos_Mbp),
        lod=max(lod)))
  }
  for(genei in genes) {
    p[[genei]] <-
      plot_feature_tbl(
        dplyr::mutate(
          dplyr::filter(gene_exon, gene==genei),
          gene = NA),
        top_snps_tbl = top_snps_tbl,
        str_rect="",
        ...) +
      ggplot2::ggtitle(genei)
  }
  if(plot_now & length(p)) {
    for(genei in genes)
      print(p[[genei]])
  }
  invisible(p)
}
#' @method autoplot gene_exon
#' @export
#' @export autoplot.gene_exon
#' @rdname plot_gene_exon
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.gene_exon <- function(x, ...)
  plot_gene_exon(x, ...)

#' @method plot gene_exon
#' @export
#' @export plot.gene_exon
#' @rdname plot_gene_exon
#' 
plot.gene_exon <- function(x, ...)
  autoplot.gene_exon(x, ...)
