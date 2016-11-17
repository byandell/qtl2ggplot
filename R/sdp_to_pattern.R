#' Convert sdp to pattern
#'
#' Convert strain distribution pattern (sdp) to letter pattern.
#'
#' @param sdp vector of sdp values
#' @param haplos letter codes for haplotypes (defaults to first 8 \code{LETTERS})
#'
#' @return vector of letter patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{sdp_to_pattern(snp)}
#'
#' @export
sdp_to_pattern <- function(sdp, haplos = LETTERS[1:8]) {
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
#' @export
#' @rdname sdp_to_pattern
sdp_to_logical <- function(sdp, haplos = LETTERS[1:8]) {
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}
