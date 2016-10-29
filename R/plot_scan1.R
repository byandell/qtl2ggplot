#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param x Output of \code{\link[qtl2scan]{scan1}}.
#'
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). Only one value allowed.
#'
#' @param chr Selected chromosomes to plot; a vector of character
#' strings.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param bgcolor Background color for the plot.
#'
#' @param altbgcolor Background color for alternate chromosomes.
#'
#' @param gap Gap between chromosomes.
#'
#' @param ... Additional graphics parameters.
#'
#' @seealso \code{\link{plot_coef}}, \code{\link{plot_coefCC}}, \code{\link{plot_snpasso}}
#'
#' @export
#' @importFrom graphics plot rect lines par axis title abline box
#'
#' @return None.
#'
#' @examples
#' # load qtl2geno package for data and genoprob calculation
#' library(qtl2geno)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, step=1, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- get_x_covar(iron)
#'
#' # perform genome scan
#' library(qtl2scan)
#' out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
#'
#' # plot the results for selected chromosomes
#' ylim <- c(0, maxlod(out)*1.02) # need to strip class to get overall max LOD
#' chr <- c(2,7,8,9,15,16)
#' plot(out, lodcolumn=1:2, chr=chr, ylim=ylim, col=c("darkslateblue","violetred"), 
#'      legend.position=c(0.1,0.9))
#'
#' # plot just one chromosome
#' plot(out, chr=8, lodcolumn=1:2, ylim=ylim, col=c("darkblue","violetred"))
#'
#' # lodcolumn can also be a column name
#' plot(out, chr=8, lodcolumn=c("liver","spleen"), ylim=ylim, col=c("darkblue","violetred"))
plot_scan1 <-
    function(x, lodcolumn=1, chr=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    # pull out map
    map <- x$map
    if(is.null(map)) stop("No map found in the input")
    if(!is.list(map)) map <- list(" "=map) # if a vector, treat it as a list with no names

    # pull out lod scores
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(x$lod))
        if(any(is.na(tmp)))
            stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(any(lodcolumn < 1 || lodcolumn > ncol(x$lod)))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(x$lod), ")")
    lod <- x$lod[,lodcolumn, drop = FALSE]

    # subset chromosomes
    if(!is.null(chr)) {
        chri <- match(chr, names(map))
        if(any(is.na(chri)))
            stop("Chromosomes ", paste(chr[is.na(chri)], collapse=", "), " not found")
        map <- map[chri]
        lod <- lod[unlist(lapply(map, names)),, drop = FALSE]
    }

    # make the plot
    ggplot_scan1(map=map, lod=lod, add=add, gap=gap,
                 bgcolor=bgcolor, altbgcolor=altbgcolor,
                 ...)
}


#' @export
#' @rdname plot_scan1
plot.scan1 <-
    function(x, lodcolumn=1, chr=NULL, add=FALSE, gap=25,
             bgcolor="gray90", altbgcolor="gray85", ...)
{
    # if snp asso result, use plot_snpasso() with just reduced snps; otherwise defaults
    if(!is.null(x$snpinfo)) {
        plot_snpasso(x, add=add, gap=gap, bgcolor=bgcolor,
                     altbgcolor=altbgcolor, ...)
    }
    else { # mostly, use plot_scan1()
        plot_scan1(x, lodcolumn=lodcolumn, chr=chr, add=add, gap=gap,
                   bgcolor=bgcolor, altbgcolor=altbgcolor, ...)
    }
}


# convert map to x-axis positions for plot_scan1
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

# boundaries of chromosomes in plot_scan1
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
