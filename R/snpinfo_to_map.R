# copied from qtl2 
snpinfo_to_map <-
  function(snpinfo)
  {
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
      stop("snpinfo$index values outside of range [1, ",
           nrow(snpinfo), "]")
    
    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)
    
    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq(along=map)) {
      u <- unique(index[[i]])
      map[[i]] <- map[[i]][u]
      names(map[[i]]) <- snp[[i]][u]
    }
    
    names(map) <- uchr
    
    map
  }
# From https://github.com/rqtl/qtl2/blob/master/R/test_util.R
is_number <- function(x) is.numeric(x) && length(x)==1
is_nonneg_number <- function(x) is_number(x) && x >= 0
is_pos_number <- function(x) is_number(x) && x > 0
