#' Summary of scan1 object
#'
#' @param object object from \code{\link[qtl2scan]{scan1}}
#' @param lodcolumn one or more lod columns
#' @param chr one or more chromosome IDs
#' @param type type of summary
#' @param drop LOD drop from maximum
#' @param show_all_snps show all SNPs if \code{TRUE}
#' @param ... other arguments not used
#'
#' @return tbl summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(summary(object))
#'
#' @export
#' @importFrom dplyr arrange desc group_by mutate summarize tbl_df ungroup
summary_scan1 <- function(object, 
                          lodcolumn=seq_len(ncol(object$lod)), 
                          chr = names(object$map), 
                          type = c("common","best"), drop=1.5, 
                          show_all_snps = TRUE,
                          ...) {
  if(is.null(object$snpinfo)) {
    # scan1 LOD summary
    chr <- as.character(chr)
    peak_pos <-
      matrix(NA, length(chr), length(lodcolumn))
    if(is.numeric(lodcolumn))
      lodcolumn <- dimnames(object$lod)[[2]][lodcolumn]
    dimnames(peak_pos) <- list(as.character(chr), lodcolumn)
    peak_lod <- peak_pos
    for(lodc in lodcolumn)
      for(chri in as.character(chr)) {
        tmp <- max(object, lodc, chri)
        if(nrow(tmp)) {
          if(is.finite(tmp[[3]])) {
            peak_pos[chri,lodc] <- mean(tmp$pos)
            ## lod column of max_scan1 has name of trait
            peak_lod[chri,lodc] <- max(tmp[[3]])
          }
        }
      }
    dplyr::tbl_df(data.frame(pheno=rep(lodcolumn, each=length(chr)),
                      chr=rep(chr, times=length(lodcolumn)),
                      pos=c(peak_pos),
                      lod=c(peak_lod)))
  } else {
    # snpinfo summary
    ## top_snps() Adapted to multiple phenotypes.
    object <- top_snps_all(object, drop, show_all_snps)
    
    if(!nrow(object))
      return(NULL)
    
    type <- match.arg(type)
    switch(type,
           best = { ## Top SNPs across all phenotypes.
             if(!nrow(object))
               return(NULL)
             dplyr::arrange(
               dplyr::mutate(object, 
                             pattern = sdp_to_pattern(sdp)),
               dplyr::desc(lod))},
           common = { ## Find most common patterns by pheno.
             dplyr::arrange(
               dplyr::mutate(
                 dplyr::ungroup(
                   dplyr::summarize(
                     dplyr::group_by(object,pheno,sdp), 
                     count=n(),
                    pct=round(100 * n() / nrow(object), 2),
                    min_lod=min(lod),
                    max_lod=max(lod),
                    max_snp=snp_id[which.max(lod)],
                    max_pos=pos_Mbp[which.max(lod)])), 
                 pattern = sdp_to_pattern(sdp)),
               dplyr::desc(max_lod))
           })
  }
}
#' @method summary scan1
#' @rdname summary_scan1
#' @export
summary.scan1 <- function(object, ...) {
  summary_scan1(object, ...)
}
