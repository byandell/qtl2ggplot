#' Top SNPs for all phenotypes scanned
#'
#' Separate fine mapping scans by allele pattern.
#'
#' @param scan1_output output of linear mixed model for \code{phename} (see \code{\link[qtl2scan]{scan1}})
#' @param drop include all SNPs within \code{drop} of max LOD (default 1.5)
#' @param show_all_snps show all SNPs if \code{TRUE}
#'
#' @return table of top_snps at maximum lod for \code{pattern}
#'
#' @export
#' @importFrom dplyr filter inner_join select
#' @importFrom tidyr gather
top_snps_all <- function (scan1_output, drop = 1.5, show_all_snps = TRUE)
{
    map <- scan1_output$map
    if (is.null(map))
        stop("No map found")
    snpinfo <- scan1_output$snpinfo
    if (is.null(snpinfo))
        stop("No snpinfo found")
    chr <- names(map)
    if (length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    ## Following is generalized from qtl2scan::top_snps()
    lod_df <- as.data.frame(subset(scan1_output, chr = chr)$lod)
    lod_df$index <- seq(nrow(lod_df))
    lod_df$snp_id <- rownames(lod_df)
    lod_df <- tidyr::gather(lod_df, pheno, lod, -snp_id, -index)
    lod_df <- dplyr::filter(lod_df, lod > max(lod, na.rm = TRUE) - drop)

    snpinfo <- snpinfo[[chr]]

    if (show_all_snps) {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -snp_id),
                                   by = "index")
    }
    else {
      snpinfo <- dplyr::inner_join(snpinfo,
                                   dplyr::select(lod_df, -index),
                                   by = "snp_id")
    }
    snpinfo
}
