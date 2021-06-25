# See https://github.com/rqtl/qtl2/blob/master/R/subset_scan1.R

# grab marker names as a vector
map2markernames <- function(map) {
  nam <- unlist(lapply(map, names))
  names(nam) <- NULL
  nam
}

# grab chromosome IDs as a vector
map2chr <- function (map) {
  chr <- rep(names(map), vapply(map, length, 0))
  names(chr) <- map2markernames(map)
  chr
}

# grab positions as a vector
map2pos <- function (map) {
  pos <- unlist(map)
  names(pos) <- map2markernames(map)
  pos
}