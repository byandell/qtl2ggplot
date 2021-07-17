# helper routine
snpinfo_to_haplos <- function(snpinfo) {
  # This routine is brittle. It depends on specific names in snpinfo and/or nc.
  alleles <- names(dplyr::select(
    snpinfo,
    -(.data$snp_id:.data$alleles)))
  # Would be better to have object that gives allele names rather than this opposite approach.
  infonames <- c("consequence","type","sdp","index","interval","on_map","pheno","lod","ensembl_gene")
  m <- match(alleles, infonames)
  alleles <- alleles[is.na(m)]
  
  # Columns in between consequence and type should be alleles.
  # If not provided, assume we are in mouse with 8.
  if((nc <- length(alleles)) < 2) {
    warning("no alleles in snpinfo; assuming 8")
    nc <- 8
  }
  LETTERS[seq_len(nc)]
}
