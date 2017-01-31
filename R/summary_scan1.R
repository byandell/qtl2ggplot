#' Summary of topsnp_pattern object
#'
#' @param object of class \code{topsnp_pattern}
#' @param type type of summary
#' @param drop LOD drop from maximum
#' @param show_all_snps show all SNPs if \code{TRUE}
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
summary_scan1 <- function(object,
                          type = c("best","common"),
                          drop = 1.5,
                          show_all_snps = FALSE,
                          ...) {
  ## Adapted to multiple phenotypes.
  object <- top_snps_all(object, drop, show_all_snps)

  if(!nrow(object))
    return(NULL)

  type <- match.arg(type)
  switch(type,
         best = { ## Top SNPs across all phenotypes.
           if(!nrow(object))
             return(NULL)
           object %>%
             mutate(pattern = sdp_to_pattern(sdp)) %>%
             arrange(desc(lod))},
         common = { ## Find most common patterns by pheno.
           object %>%
             group_by(pheno,sdp) %>%
             summarize(count=n(),
                       pct=round(100 * n() / nrow(object), 2),
                       min_lod=min(lod),
                       max_lod=max(lod),
                       max_snp=snp_id[which.max(lod)],
                       max_pos=pos_Mbp[which.max(lod)]) %>%
             ungroup() %>%
             mutate(pattern = sdp_to_pattern(sdp)) %>%
             arrange(desc(max_lod))
         })
}
#' @method summary scan1
#' @rdname summary_scan1
#' @export
summary.scan1 <- function(object, type, ...) {
  summary_scan1(object, type, ...)
}
