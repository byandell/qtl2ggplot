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
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). One or more value(s) allowed.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param drop.hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#'
#' @param col.hilit Color of highlighted points
#'
#' @param col Color of other points
#'
#' @param gap Gap between chromosomes.
#'
#' @param ylim y-axis limits
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}.
#' @param legend.position Position of legend (default \code{"none"}).
#'
#' @param ... Additional graphics parameters.
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
#' # download snp info from web
#' tmpfile <- tempfile()
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/c2_snpinfo.rds")
#' download.file(file, tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#'
#' # calculate strain distribution patterns
#' snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
#'
#' # switch map in allele probabilities to Mbp
#' apr$map <- DOex$pmap
#'
#' # identify equivalent snps
#' snpinfo <- index_snps(DOex$pmap, snpinfo)
#'
#' # convert to snp probabilities
#' snppr <- genoprob_to_snpprob(apr, snpinfo)
#'
#' # perform SNP association analysis (here, ignoring residual kinship)
#' out_snps <- scan1(snppr, DOex$pheno)
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
#' autoplot(out_snps, snpinfo, drop.hilit=1.5)
#'
#' # highlight SDP patterns in SNPs; connect with lines.
#' autoplot(out_snps, snpinfo, patterns="all",drop.hilit=4)
#'
#' # highlight top SDP patterns in SNPs; connect with lines.
#' autoplot(out_snps, snpinfo, patterns="hilit",drop.hilit=4)
#' }
#'
#' @seealso \code{\link{ggplot_scan1}}, \code{\link{ggplot_coef}}, \code{\link{ggplot_coefCC}}
#'
#' @export
#'
plot_snpasso <-
    function(scan1output, snpinfo, lodcolumn=1, show_all_snps=TRUE, drop.hilit=NA,
             col.hilit="violetred", col="darkslateblue",
             ylim=NULL, gap=25,
             bgcolor="gray90", altbgcolor="gray85",
             ...)
{
    uindex <- unique(snpinfo$index)
    if(length(uindex) != nrow(scan1output))
        stop("Something is wrong with snpinfo$index.\n",
             "      length(unique(snpinfo$index)) [",
             length(unique(snpinfo$index)), "] != nrow(scan1output) [",
             nrow(scan1output), "].")

    if(any(snpinfo$index[uindex] != uindex))
        stop("Something is wrong with snpinfo$index.\n",
             "      snpinfo$index[u] should == u for values in snpinfo$index")

    plot_snpasso_internal(scan1output, snpinfo, lodcolumn, show_all_snps, drop.hilit,
                          col.hilit, col,
                          ylim, gap,
                          bgcolor, altbgcolor,
                          ...)
}

plot_snpasso_internal <- function(scan1output, snpinfo, lodcolumn, show_all_snps, drop.hilit,
                                  col.hilit, col,
                                  ylim, gap,
                                  bgcolor, altbgcolor,
                                  patterns = c("none","all","hilit"),
                                  lines = (patterns != "none"), points = TRUE,
                                  legend.position = ifelse((patterns != "none"), "right", "none"),
                                  legend.title = ifelse((patterns != "none"), "pattern", "pheno"),
                                  reorder = TRUE,
                                  ...) {

  map <- snpinfo_to_map(snpinfo)
  
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
  if(patterns != "none")
    pattern <- sdp_to_pattern(snpinfo$sdp)

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
                                  col.hilit, drop.hilit, maxlod)
  # settings$pattern will be either SDP patterns or thresholding by drop.hilit.

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
# taken from CCSanger::sdp_to_pattern
sdp_to_pattern <- function(sdp, haplos = LETTERS[1:8]) {
  sapply(sdp, function(x, haplos) {
    ref <- as.logical(intToBits(x)[seq_along(haplos)])
    paste(paste(haplos[!ref], collapse = ""),
          paste(haplos[ref], collapse = ""),
          sep = ":")
  }, haplos)
}
