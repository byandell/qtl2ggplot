#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param x Output of \code{\link[qtl2]{scan1}}.
#'
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2]{insert_pseudomarkers}}.
#'
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). One or more value(s) allowed.
#'
#' @param chr Selected chromosomes to plot; a vector of character
#' strings.
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param gap Gap between chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @seealso \code{\link{ggplot_coef}}, \code{\link{ggplot_snpasso}}
#'
#' @export
#' @importFrom graphics plot rect lines par axis title abline box
#' @importFrom qtl2 align_scan1_map subset_scan1
#'
#' @return None.
#'
#' @examples
#' # load qtl2 package for data and genoprob calculation
#' library(qtl2)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # plot the results for selected chromosomes
#' chr <- c(2,7,8,9,15,16)
#' ggplot_scan1(out, map, lodcolumn=1:2, chr=chr, col=c("darkslateblue","violetred"),
#'      legend.position=c(0.1,0.9))
#'
#' # plot just one chromosome
#' ggplot_scan1(out, map, chr=8, lodcolumn=1:2, col=c("darkblue","violetred"))
#'
#' # can also use autoplot from ggplot2
#' # lodcolumn can also be a column name
#' library(ggplot2)
#' autoplot(out, map, chr=8, lodcolumn=c("liver","spleen"), col=c("darkblue","violetred"))
ggplot_scan1 <-
    function(x, map, lodcolumn=1, chr=NULL, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    if(!is.list(map)) map <- list(" " = map) # if a vector, treat it as a list with no names
    
    # subset chromosomes
    if(!is.null(chr)) {
      chri <- match(chr, names(map))
      if(any(is.na(chri)))
        stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
      x <- qtl2::subset_scan1(x, map, chr)
      map <- map[chri]
    }
    
    # align scan1 output and map
    tmp <- qtl2::align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map
    
    if(nrow(x) != length(unlist(map)))
        stop("nrow(x) [", nrow(x), "] != number of positions in map [",
             length(unlist(map)), "]")

    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(x))
        if(any(is.na(tmp)))
            stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(any(lodcolumn < 1) || any(lodcolumn > ncol(x)))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x), ")")
    lod <- unclass(x)[,lodcolumn, drop = FALSE]

    # subset chromosomes
    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        map <- map[chri]
        lod <- lod[unlist(lapply(map, names)),, drop = FALSE]
    }

    # make the plot
    ggplot_scan1_internal(map=map, lod=lod, gap=gap,
                 bgcolor=bgcolor, altbgcolor=altbgcolor,
                 ...)
}




# convert map to x-axis positions for ggplot_scan1
map_to_xpos <-
    function(map, gap)
{
    if(length(map)==1) return(map[[1]])

    chr_range <- vapply(map, range, c(0,1), na.rm=TRUE)

    result <- map[[1]]-chr_range[1,1] + gap/2
    for(i in 2:length(map)) {
        result <- c(result,
                    map[[i]] - chr_range[1,i] + gap + max(result, na.rm=TRUE))
    }
    result
}

# boundaries of chromosomes in ggplot_scan1
# first row: left edges
# second row: right edges
map_to_boundaries <-
    function(map, gap)
{
    if(length(map)==1)
        return(cbind(range(map[[1]], na.rm=TRUE)))

    # range of each chromosome
    chr_range <- lapply(map, range, na.rm=TRUE)

    # corresponding xpos, as matrix with two rows
    startend <- matrix(map_to_xpos(chr_range, gap), nrow=2)

    startend[1,] <- startend[1,] - gap/2
    startend[2,] <- startend[2,] + gap/2

    startend
}

#' @export autoplot.scan1
#' @export
#' @method autoplot scan1
#' @rdname ggplot_scan1
#'
#' @importFrom ggplot2 autoplot
#'
autoplot.scan1 <-
  function(x, map, lodcolumn=1, chr=NULL, gap=25,
           bgcolor="gray90", altbgcolor="gray85", ...)
  {
    # if snp asso result, use ggplot_snpasso() with just reduced snps; otherwise defaults
    if(is.data.frame(map) && "index" %in% names(map)) {
      ggplot_snpasso(x, snpinfo=map, lodcolumn=lodcolumn, gap=gap, bgcolor=bgcolor,
                   altbgcolor=altbgcolor, ...)
    }
    else { # mostly, use ggplot_scan1()
      ggplot_scan1(x, map=map, lodcolumn=lodcolumn, chr=chr, gap=gap,
                 bgcolor=bgcolor, altbgcolor=altbgcolor, ...)
    }
  }

# convert map to list of indexes to LOD vector
map_to_index <-
    function(map)
{
    if(length(map)==1) {
        map[[1]] <- seq(along=map[[1]])
        return(map)
    }

    lengths <- vapply(map, length, 0)
    split(1:sum(lengths), rep(seq(along=map), lengths))
}
