#' Plot gene locations for a genomic interval
#'
#' Plot gene locations for a genomic interval, as rectangles with gene
#' symbol (and arrow indicating strand/direction) below.
#'
#' @param genes Data frame containing \code{start} and \code{stop} in
#' bp, \code{strand} (as \code{"-"}, \code{"+"}, or \code{NA}), and
#' \code{Name}.
#' @param xlim x-axis limits (in Mbp)
#' @param minrow Minimum number of rows of genes
#' @param padding Proportion to pad with white space around the genes
#' @param colors Vectors of colors, used sequentially and then re-used.
#' @param x Object of class \code{genes}
#' @param ... Optional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return None.
#'
#' @keywords hgraphics
#' @export
#' @importFrom graphics strheight strwidth text plot par rect abline box
#'
#' @examples
#' genes <- data.frame(chr = c("6", "6", "6", "6", "6", "6", "6", "6"),
#'                     start = c(139.988753, 140.680185, 141.708118, 142.234227, 142.587862,
#'                               143.232344, 144.398099, 144.993835),
#'                     stop  = c(140.041457, 140.826797, 141.773810, 142.322981, 142.702315,
#'                               143.260627, 144.399821, 145.076184),
#'                     strand = c("-", "+", "-", "-", "-", NA, "+", "-"),
#'                     Name = c("Plcz1", "Gm30215", "Gm5724", "Slco1a5", "Abcc9",
#'                              "4930407I02Rik", "Gm31777", "Bcat1"),
#'                     stringsAsFactors=FALSE)
#' ggplot_genes(genes, xlim=c(140, 146))

# create an empty plot with test x- and y-axis limits
ggplot_genes <-
    function(genes, xlim=NULL, minrow=4, padding=0.2,
             colors=c("black", "red3", "green4", "blue3", "orange"),
             ...)
{
    # need both 'start' and 'stop' columns with no missing values
    stopifnot(!any(is.na(genes$start)), !any(is.na(genes$stop)))

    # make sure genes are ordered by their start values
    if(any(diff(genes$start) < 0))
        genes <- genes[order(genes$start, genes$stop),]

    # grab data
    start <- genes$start # assume in Mbp
    end <- genes$stop    # assume in Mbp
    strand <- as.character(genes$strand)
    name <- as.character(genes$Name)

    # missing names: use ?
    name[is.na(name)] <- "?"
    
    if(set_xlim <- is.null(xlim)) {
      xlim <- range(c(start, end), na.rm=TRUE)
    } else {
      # Keep only if start and end are in range. 
      keep <- end >= xlim[1] & start <= xlim[2]
      if(!any(keep)) {
        cat("no genes within xlim interval\n", file = stderr())
        return(NULL)
      }
      start <- start[keep]
      end <- end[keep]
      strand <- strand[keep]
      name <- name[keep]
    }
    
    # arrow annotation re direction, to place after gene name
    dir_symbol <- rep('', length(name))
    right <- !is.na(strand) & strand == "+"
    if(any(right))
      dir_symbol[right] <- "~phantom('') %->% phantom('')"
    left <- !is.na(strand) & strand == "-"
    if(any(left))
      dir_symbol[left] <- "~phantom('') %<-% phantom('')"
    
    # initial determination of text size
    text_cex <- 1
    maxy <- minrow
    height <- 1/maxy
    
    # ggplot2 does not allocate space for text, so this is approximate
    text_size <- 4 # default text size
    plot_width <- 6 # default plot width
    str_width <- function(chars, text_fudge = 4) {
      strwidth(chars, "inches") * diff(xlim) * text_size /
        (text_fudge * plot_width)
    }
    
    # horizontal padding
    space <- str_width(' ')
    # end of strings
    end_str <- end + space + str_width(name) + 
      str_width(expression(dir_symbol))
    # adjust text size and determine vertical location of genes
    for(it in 1:2) { # go through all of this twice
        # figure out how to arrange genes vertically
        #     + number of rows of genes
        # (function defined in src/arrange_genes.cpp)
        y <- arrange_genes(start, end_str)

        maxy <- max(c(y, minrow))
        height <- 1/maxy
    }
    if(set_xlim) {
      # make room for names
      xlim[1] <- xlim[1] - str_width('  ')
      xlim[2] <- max(end_str)
    }

    ypos <- seq(height/2, by=height, length=maxy)
    y <- ypos[y]
    rect_height <- height*(1-padding)
    rect_top <- y - rect_height/2
    rect_bottom <- y + rect_height/2

    colors <- rep(colors, length = length(y))
    
    ggplot_genes_internal(start, end, strand, rect_top, rect_bottom, 
                 colors, space, y, dir_symbol, name, xlim, ...)
}

#' @method autoplot genes
#' @export
#' @export autoplot.genes
#' @rdname ggplot_genes
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.genes <- function(x, ...)
  ggplot_genes(x, ...)
