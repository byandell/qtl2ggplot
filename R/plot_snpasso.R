#' Plot SNP associations
#'
#' Plot SNP associations, with possible expansion from distinct snps to all snps.
#'
#' @param scan1output Output of \code{\link[qtl2scan]{scan1}}.  Should
#' contain an attribute, \code{"snpinfo"}, as when
#' \code{\link[qtl2scan]{scan1}} are run with SNP probabilities
#' produced by \code{\link[qtl2scan]{genoprob_to_snpprob}}.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param cex Character expansion for the points (default 0.5)
#'
#' @param pch Plotting character for the points (default 16)
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
#' library(qtl2geno)
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
#' library(qtl2scan)
#' snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
#'
#' # switch map in allele probabilities to Mbp
#' apr$map <- DOex$pmap
#'
#' # convert to snp probabilities
#' snppr <- genoprob_to_snpprob(apr, snpinfo)
#'
#' # perform SNP association analysis (here, ignoring residual kinship)
#' out_snps <- scan1(snppr, DOex$pheno)
#'
#' # plot results
#' library(qtl2ggplot)
#' plot_snpasso(out_snps)
#'
#' # can also just type plot()
#' plot(out_snps)
#'
#' # plot just subset of distinct SNPs
#' plot_snpasso(out_snps, show_all_snps=FALSE)
#'
#' # highlight the top snps (with LOD within 1.5 of max)
#' plot(out_snps, drop.hilit=1.5)
#' 
#' # highlight SDP patterns in SNPs; connect with lines.
#' plot(out_snps, patterns="all",drop.hilit=4,cex=2)
#' 
#' # highlight top SDP patterns in SNPs; connect with lines.
#' plot(out_snps, patterns="hilit",drop.hilit=4,cex=2)
#' }
#'
#' @seealso \code{\link{plot_scan1}}, \code{\link{plot_coef}}, \code{\link{plot_coefCC}}
#' 
#' @export
#'
plot_snpasso <-
    function(scan1output, show_all_snps=TRUE, drop.hilit=NA,
             col.hilit="violetred", col="darkslateblue",
             pch=16, cex=0.5, ylim=NULL, gap=25,
             bgcolor="gray90", altbgcolor="gray85",
             ...)
{
    plot_snpasso_internal(scan1output, show_all_snps, drop.hilit,
                          col.hilit, col,
                          pch, cex, ylim, gap,
                          bgcolor, altbgcolor,
                          ...)
}

plot_snpasso_internal <- function(scan1output, show_all_snps, drop.hilit,
                                  col.hilit, col,
                                  pch, cex, ylim, gap,
                                  bgcolor, altbgcolor,
                                  patterns = c("none","all","hilit"),
                                  lines = (patterns != "none"),
                                  legend.position = ifelse((patterns != "none"), "right", "none"), 
                                  legend.title = ifelse((patterns != "none"), "pattern", "pheno"),
                                  ...) {
  patterns <- match.arg(patterns)
  if(patterns != "none")
    group <- sdp_to_pattern(scan1output$snpinfo[[1]]$sdp)
  
  if(show_all_snps)
    scan1output <- expand_snp_results(scan1output)
  
  # maximum LOD
  maxlod <- max(scan1output$lod[,1], na.rm=TRUE)
  
  if(is.null(ylim))
    ylim <- c(0, maxlod*1.02)
  
  settings <- color_patterns_set(scan1output, patterns,
                                 col, group, show_all_snps, 
                                 col.hilit, drop.hilit, maxlod)
  
  plot_scan1(scan1output, lodcolumn=1, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
             gap=gap, cex=cex, pch=pch, 
             col = settings$col,
             group = settings$group,
             lines = lines, points = TRUE, 
             legend.position = legend.position,
             legend.title = legend.title,
             patterns = patterns,
             ...)
}


# expand snp association results according to snpinfo
expand_snp_results <-
    function(snp_results)
{
    snpinfo <- snp_results$snpinfo
    if(is.null(snpinfo)) stop("No snpinfo found")
    map <- snp_results$map
    if(is.null(map)) stop("No map found")

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(length(snp_results$lod) != length(unlist(map)))
        stop("length(snp_results$lod) [", length(snp_results$lod), "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    lodindex <- split(seq(along=snp_results$lod), rep(names(map), vapply(map, length, 0)))

    result <- NULL
    for(i in seq(along=map)) {
        map[[i]] <- snpinfo[[i]]$pos
        names(map[[i]]) <- snpinfo[[i]]$snp
        result <- rbind(result, snp_results$lod[lodindex[[i]],,drop=FALSE][snpinfo[[i]]$index,,drop=FALSE])
    }

    list(lod=result,
         map=map)
}
