#' Convert sdp to pattern
#'
#' Convert strain distribution pattern (sdp) to letter pattern.
#' Taken from package `qtl2pattern` for internal use here.
#'
#' @param sdp vector of sdp values
#' @param haplos letter codes for haplotypes (required)
#'
#' @return vector of letter patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @importFrom assertthat assert_that
#' 
sdp_to_pattern <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
