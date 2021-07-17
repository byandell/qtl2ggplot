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
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Download SNP info for DOex from web and read as RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Convert to SNP probabilities
#' snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
#' snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
#' 
#' # Scan SNPs.
#' scan_snppr <- qtl2::scan1(snppr, DOex$pheno)
#'
#' # plot results
#' ggplot_snpasso(scan_snppr, snpinfo, drop_hilit=1.5)
#'
#' # can also just type autoplot() if ggplot2 attached
#' library(ggplot2)
#'
#' # plot just subset of distinct SNPs
#' autoplot(scan_snppr, snpinfo, show_all_snps=FALSE, drop_hilit=1.5)
#'
#' # highlight SDP patterns in SNPs; connect with lines.
#' autoplot(scan_snppr, snpinfo, patterns="all",drop_hilit=4)
#'
#' # query function for finding genes in region
#' gene_dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
#' query_genes <- qtl2::create_gene_query_func(gene_dbfile)
#' genes <- query_genes(2, 97, 98)
#'
#' # plot SNP association results with gene locations
#' autoplot(scan_snppr, snpinfo, patterns="hilit", drop_hilit=1.5, genes=genes)
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
                                  legend.title = ifelse((patterns != "none"), "pattern", "pheno"),
                                  reorder = TRUE,
                                  haplos = snpinfo_to_haplos(snpinfo),
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

