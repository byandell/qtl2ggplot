# helper routine
sdp_to_pattern <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
