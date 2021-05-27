#' Plot SNP associations
#'
#' Plot SNP associations, with possible expansion from distinct snps to all snps.
#'
#' @param scan1output Output of \code{\link[qtl2]{scan1}}.  Should
#' contain an attribute, \code{"snpinfo"}, as when
#' \code{\link[qtl2]{scan1}} are run with SNP probabilities
#' produced by \code{\link[qtl2]{genoprob_to_snpprob}}.
#'
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived from with
#'     \code{\link[qtl2]{index_snps}}):
#' \itemize{
#' \item \code{chr} - Character string or factor with chromosome
#' \item \code{pos} - Position (in same units as in the \code{"map"}
#'     attribute in \code{genoprobs}.
#' \item \code{sdp} - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' \item \code{snp} - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' \item \code{index} - Indices that indicate equivalent
#'     groups of SNPs.
#' \item \code{intervals} - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' \item \code{on_map} - Indicate whether SNP coincides with a marker
#'     in the \code{genoprobs}
#' }
#'
#' @param genes Optional data frame containing gene information for
#' the region, with columns `start` and `stop` in Mbp, `strand`
#' (as `"-"`, `"+"`, or `NA`), and `Name`. If included, a
#' two-panel plot is produced, with SNP associations above and
#' gene locations below.
#'
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). One or more value(s) allowed.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param drop_hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#'
#' @param col_hilit Color of highlighted points
#'
#' @param col Color of other points
#'
#' @param ylim y-axis limits
#'
#' @param gap Gap between chromosomes.
#'
#' @param minlod Minimum LOD to display. (Mostly for GWAS, in which
#'     case using `minlod=1` will greatly increase the plotting speed,
#'     since the vast majority of points would be omittted.
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#' 
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color and `altbgcolor`
#' to control the background color on alternate chromosomes.
#' `cex` for character expansion for the points (default 0.5),
#' `pch` for the plotting character for the points (default 16),
#' and `ylim` for y-axis limits.
#'
#' @examples
#' \dontrun{
#' # load example DO data from web
#' library(qtl2)
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#'
#' # subset to chr 2
#' DOex <- DOex[,"2"]
#'
#' # calculate genotype probabilities and convert to allele probabilities
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' apr <- genoprob_to_alleleprob(pr)
#'
#' # query function for grabbing info about variants in region
#' snp_dbfile <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
#' query_variants <- create_variant_query_func(snp_dbfile)
#'
#' # SNP association scan
#' out_snps <- scan1snps(apr, DOex$pmap, DOex$pheno, query_func=query_variants,
#'                       chr=2, start=97, end=98, keep_all_snps=TRUE)
#'
#' # plot results
#' ggplot_snpasso(out_snps, snpinfo)
#'
#' # can also just type autoplot() if ggplot2 attached
#' library(ggplot2)
#' autoplot(out_snps, snpinfo)
#'
#' # plot just subset of distinct SNPs
#' autoplot(out_snps, snpinfo, show_all_snps=FALSE)
#'
#' # highlight the top snps (with LOD within 1.5 of max)
#' autoplot(out_snps, snpinfo, drop_hilit=1.5)
#'
#' # highlight SDP patterns in SNPs; connect with lines.
#' autoplot(out_snps, snpinfo, patterns="all",drop_hilit=4)
#'
#' # highlight top SDP patterns in SNPs; connect with lines.
#' autoplot(out_snps, snpinfo, patterns="hilit",drop_hilit=4)
#'
#' # query function for finding genes in region
#' gene_dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
#' query_genes <- create_gene_query_func(gene_dbfile)
#' genes <- query_genes(2, 97, 98)
#'
#' # plot SNP association results with gene locations
#' autoplot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)
#' }
#'
#' @seealso \code{\link{ggplot_scan1}}, \code{\link{ggplot_coef}}
#'
#' @export
#' @importFrom assertthat assert_that
#'
ggplot_snpasso <-
    function(scan1output, snpinfo, genes=NULL, lodcolumn=1, show_all_snps=TRUE, drop_hilit=NA,
             col_hilit="violetred", col="darkslateblue",
             ylim=NULL, gap=25, minlod=0,
             bgcolor="gray90", altbgcolor="gray85",
             ...)
{
    if(is.null(scan1output)) stop("scan1output is NULL")
    if(is.null(snpinfo)) stop("snpinfo is NULL")
    
    if(!is_nonneg_number(gap)) stop("gap should be a single non-negative number")
    if(!is_nonneg_number(minlod)) stop("minlod should be a single non-negative number")
    
    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
#    if(length(lodcolumn) > 1) { # If length > 1, take first value
#      warning("lodcolumn should have length 1; only first element used.")
#      lodcolumn <- lodcolumn[1]
#    }
    if(is.character(lodcolumn)) { # turn column name into integer
      tmp <- match(lodcolumn, colnames(scan1output))
      if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
      lodcolumn <- tmp
    }
    if(any(lodcolumn < 1) || any(lodcolumn > ncol(scan1output)))
      stop("lodcolumn [", paste(lodcolumn, collapse = ","),
           "] out of range (should be in 1, ..., ", ncol(scan1output), ")")
    scan1output <- scan1output[,lodcolumn,drop=FALSE]
    
    if(!is.null(genes)) {
      if(length(unique(snpinfo$chr)) > 1) {
        warning("genes ignored if there are multiple chromosomes")
      } else {
        return(ggplot_snpasso_and_genes(scan1output, snpinfo, show_all_snps=show_all_snps,
                                       drop_hilit=drop_hilit, col_hilit=col_hilit,
                                       col=col, gap=gap, minlod=minlod, genes=genes, ...) )
      }
    }
    
    if(nrow(scan1output) == nrow(snpinfo) && all(rownames(scan1output) == snpinfo$snp)) {
      show_all_snps <- FALSE
      
      # Keep rows with at least one LOD above minlod.
      keep <- apply(scan1output, 1, function(x) any(x>=minlod))
      snpinfo <- snpinfo[keep,, drop=FALSE]
      scan1output <- scan1output[keep,, drop=FALSE]
    }
    else {
      snpinfo_spl <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))
      
      uindex <- unlist(lapply(snpinfo_spl, function(a) unique(a$index)))
      if(length(uindex) != nrow(scan1output)) {
        stop("Something is wrong with snpinfo$index.\n",
             "      no. unique indexes [",
             length(uindex), "] != nrow(scan1output) [",
             nrow(scan1output), "].")
      }
      
      for(i in seq_along(snpinfo_spl)) {
        uindex <- unique(snpinfo_spl[[i]]$index)
        if(any(snpinfo_spl[[i]]$index[uindex] != uindex)) {
          stop("Something is wrong with snpinfo$index.\n",
               "      on each chr, index[u] should == u for each index u")
        }
      }
    }
    
    # set values < minlod to NA so they're not plotted
    scan1output[scan1output < minlod] <- NA
    
    # internal function to give defaults to hidden graphics parameters
    ggplot_snpasso_internal(scan1output, snpinfo, lodcolumn, show_all_snps, drop_hilit,
                          col_hilit, col,
                          ylim, gap,
                          bgcolor, altbgcolor,
                          ...)
}

ggplot_snpasso_internal <- function(scan1output, snpinfo, lodcolumn, show_all_snps, drop_hilit,
                                  col_hilit, col,
                                  ylim, gap,
                                  bgcolor, altbgcolor,
                                  patterns = c("none","all","hilit"),
                                  lines = (patterns != "none"), points = TRUE,
                                  legend.position = ifelse((patterns != "none"), "right", "none"),
                                  legend.title = ifelse((patterns != "none"), "pattern", "pheno"),
                                  reorder = TRUE,
                                  ...) {

  # subset on lodcolumns
  scan1output <- subset(scan1output, lodcolumn = lodcolumn)
  lodcolumn <- seq_len(ncol(scan1output))
  
  # reorder columns of scan1output by decreasing LOD
  if(reorder) {
    o <- order(-apply(scan1output, 2, max))
    scan1output <- modify_object(scan1output,
                                 scan1output[, o, drop=FALSE])
  }

  patterns <- match.arg(patterns)
  if(patterns != "none") {
    haplos <- snpinfo_to_haplos(snpinfo)
    pattern <- sdp_to_pattern(snpinfo$sdp, haplos)
  }

  map <- snpinfo_to_map(snpinfo)

  if(show_all_snps) {
    tmp <- expand_snp_results(scan1output, map, snpinfo)
    scan1output <- modify_object(scan1output, tmp$lod)
    map <- tmp$map
  }
  
  # maximum LOD
  maxlod <- max(unclass(scan1output), na.rm=TRUE)

  if(is.null(ylim))
    ylim <- c(max(0, min(unclass(scan1output), na.rm=TRUE)),
              maxlod*1.02)

  settings <- color_patterns_set(scan1output, snpinfo, patterns,
                                  col, pattern, show_all_snps,
                                  col_hilit, drop_hilit, maxlod)
  # settings$pattern will be either SDP patterns or thresholding by drop_hilit.

  ggplot_scan1(scan1output, map=map, lodcolumn=lodcolumn, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
             gap = gap,
             col = settings$col,
             pattern = settings$pattern,
             shape = settings$shape,
             lines = lines, points = points,
             legend.position = legend.position,
             legend.title = legend.title,
             patterns = patterns,
             ...)
}


# expand snp association results according to snpinfo
expand_snp_results <-
    function(snp_results, map, snpinfo)
{
    snpinfo <- split(snpinfo, snpinfo$chr)

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(nrow(snp_results) != length(unlist(map)))
        stop("nrow(snp_results) [", nrow(snp_results), "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    lodindex <- split(1:nrow(snp_results), rep(names(map), vapply(map, length, 0)))

    result <- NULL
    for(i in seq(along=map)) {
        revindex <- rev_snp_index(snpinfo[[i]])

        map[[i]] <- snpinfo[[i]]$pos
        names(map[[i]]) <- snpinfo[[i]]$snp
        result <- rbind(result, unclass(snp_results)[lodindex[[i]],,drop=FALSE][revindex,,drop=FALSE])
        rownames(result) <- snpinfo[[i]]$snp
    }

    list(lod=result,
         map=map)
}

# reverse index
rev_snp_index <-
    function(snpinfo)
{
    index_spl <- split(1:nrow(snpinfo), snpinfo$index)
    revindex <- rep(seq(along=index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex

    revindex
}

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
# taken from qtl2pattern::sdp_to_pattern
sdp_to_pattern <- function(sdp, haplos) {
  assertthat::assert_that(!missing(haplos))
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
# taken from qtl2pattern:::snpinfo_to_haplos
snpinfo_to_haplos <- function(snpinfo) {
  # This routine is brittle. It depends on specific names in snpinfo and/or nc.
  snp_id <- NULL # trick R check.
  alleles <- names(dplyr::select(
    snpinfo,
    -(snp_id:alleles)))
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
# From https://github.com/rqtl/qtl2/blob/master/R/test_util.R
is_number <- function(x) is.numeric(x) && length(x)==1
is_nonneg_number <- function(x) is_number(x) && x >= 0
is_pos_number <- function(x) is_number(x) && x > 0
