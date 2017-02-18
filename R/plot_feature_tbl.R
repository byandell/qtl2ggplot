#' Staggered plot of genes in interval
#'
#' Plot genes as rectangles followed by names.
#' Stagger genes for easy reading.
#' Written original by Dan Gatti 2013-02-13
#'
#' @param x tbl of gene information from \code{\link{get_mgi_features}}
#' @param rect_col fill color of rectangle (default "grey70")
#' @param strand_col edge color of rectangle by strand from \code{x} (default -="blue", +="red"; none if NULL)
#' @param type_col color of type from \code{x} (default "black" for gene, "blue" for pseudogene; none if NULL)
#' @param text_size size of text (default 3)
#' @param xlim horizontal axis limits (default is range of features)
#' @param snp_pos position of SNPs in bp if used (default NULL)
#' @param snp_lod LOD of SNPs (for color plotting)
#' @param top_snps_tbl table from \code{\link[qtl2scan]{top_snps}}
#' @param snp_col color of SNP vertical lines (default "grey70")
#' @param extend extend region for SNPs in bp (default 5000)
#' @param ... additional arguments (not used)
#'
#' @return data frame of gene information (invisible)
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#'    Daniel Gatti, \email{Dan.Gatti@@jax.org}
#' @references \url{https://github.com/dmgatti/DOQTL/blob/master/R/gene.plot.R}
#' @keywords hplot
#'
#' @export
#' 
#' @importFrom ggplot2 aes element_blank geom_rect geom_text geom_vline
#' ggplot scale_color_gradient theme xlab ylab
#' @importFrom dplyr arrange filter
#' @importFrom CCSanger convert_bp
#' 
plot_feature_tbl <- function(x,
                             rect_col = "grey70",
                             strand_col = c("-"="#1b9e77", "+"="#d95f02"),
                             type_col = c(gene="black", pseudogene="#1b9e77", other="#d95f02"),
                             text_size = 3,
                             xlim = NULL,
                             snp_pos = top_snps_tbl$pos_Mbp,
                             snp_lod = top_snps_tbl$lod,
                             top_snps_tbl = NULL,
                             snp_col = "grey70",
                             extend = 5000,
                             ...) {
  # If we have no genes, just plot an empty frame and return.
  if(is.null(x) || length(x) == 0) {
    plot(0, 0, col = 0, xlab = "", xaxs = "i", ylab = "", yaxt = "n", ...)
    return()
  } # if(is.null(x) || nrow(x) == 0)
  
  x$start <- CCSanger::convert_bp(x$start)
  x$stop <- CCSanger::convert_bp(x$stop)
  x <- dplyr::arrange(x, desc(type),strand,start)
  
  # Expand rect_col and text_col; add Name.
  if(is.null(rect_col))
    rect_col <- "grey70"
  rect_col <- rep_len(rect_col, nrow(x))
  if(is.null(strand_col)) {
    rect_edge <- "grey30"
  } else {
    rect_edge <- strand_col[match(x$strand,
                                  names(strand_col))]
    ## Fill in missing values with fill color.
    miss_edge <- is.na(rect_edge)
    if(any(miss_edge))
      rect_edge[miss_edge] <- "grey30"
  }
  if(is.null(type_col)) {
    text_size <- 0
  } else {
    text_col <- type_col[match(x$type,
                               names(type_col),
                               nomatch=3)]
  }
  
  # Subset data to plot limits.
  if(is.null(xlim))
    xlim <- c(min(x$start), max(x$stop))
  else {
    x <- dplyr::filter(x, stop >= xlim[1] & start <= xlim[2])
    # If we have no genes to plot, just return.
    if(nrow(x) == 0) {
      warning("no genes in interval")
      return()
    }
  }
  retval <- get.gene.locations(x, xlim, text_size, ...)
  
  x$Name[is.na(x$Name)] <- ""
  
  # Plot the genes.
  x$bottom <- -retval$bottom
  nudge <- retval$nudge
  ## Offset between rectangles.
  rowheight <- 1
  boxheight = rowheight / 1.15
  offset = 0.15 * boxheight
  
  ## Arrange snp_pos by in decreasing snp_lod
  if(!is.null(snp_pos) & !is.null(snp_lod)) {
    o <- order(snp_lod)
    snp_pos <- snp_pos[o]
    snp_lod <- snp_lod[o]
  }
  
  p <- ggplot2::ggplot(x) +
    #    scale_y_continuous(expand = c(0,0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())
  snp_vline <- function(p, snp_pos, snp_lod, keep, xlim, extend) {
    if(is.numeric(snp_pos) & !is.null(snp_lod)) {
      
      snp_pos <- CCSanger::convert_bp(snp_pos)
      keep <- snp_pos >= xlim[1] - extend &
        snp_pos <= xlim[2] + extend
      if(any(keep)) {
        if(is.null(snp_lod)) {
          snp_col <- rep(snp_col, length = length(snp_pos))
          p + ggplot2::geom_vline(xintercept = CCSanger::convert_bp(snp_pos[keep], FALSE),
                                  linetype = "dashed", col = snp_col[keep])
        } else {
          snp_col <- snp_lod
          tmp <- data.frame(pos = CCSanger::convert_bp(snp_pos[keep], FALSE),
                            lod = snp_col[keep])
          p + ggplot2::geom_vline(data = tmp,
                                  ggplot2::aes(xintercept = pos, col = lod),
                                  linetype = "dashed") +
            ggplot2::scale_color_gradient(low = "grey90", high = "grey10")
        }
      }
    } else
      p
  }
  if(!is.null(snp_pos) & nrow(x) > 1) {
    p <- snp_vline(p, snp_pos, snp_lod, keep, xlim, extend)
  }
  p <- p +
    ggplot2::geom_rect(
      mapping = ggplot2::aes(xmin = CCSanger::convert_bp(start, FALSE),
                             xmax = CCSanger::convert_bp(stop, FALSE),
                             ymin = bottom,
                             ymax = bottom - 1 + offset),
      fill = rect_col,
      color = rect_edge) +
    ggplot2::xlab(paste("Chr", x$seqid[1], "(Mb)")) +
    ggplot2::ylab("")
  if(!is.null(snp_pos) & nrow(x) == 1) {
    ## If only one Gene, then put SNP dashes in front of rectangles to show overlap.
    p <- snp_vline(p, snp_pos, snp_lod, keep, xlim, extend)
  }
  if(!is.null(type_col)) {
    ## Want to remove entries with no Name.
    ## Do as new data here?
    if(!all(x$Name=="")) {
      p <- p +
        ggplot2::geom_text(
          mapping = ggplot2::aes(x = CCSanger::convert_bp((stop + nudge), FALSE),
                                 y = bottom - 0.5 + offset,
                                 label = Name,
                                 hjust = 0,
                                 vjust = 0.5),
          size = text_size,
          color = text_col)
    }
  }
  p
}
#' @method autoplot feature_tbl
#' @export
#' @export autoplot.feature_tbl
#' @rdname plot_feature_tbl
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.feature_tbl <- function(x, ...)
  plot_feature_tbl(x, ...)

#' @method plot feature_tbl
#' @export
#' @export plot.feature_tbl
#' @rdname plot_feature_tbl
#' 
plot.feature_tbl <- function(x, ...)
  autoplot.feature_tbl(x, ...)

#' Helper function to set gene locations on plot.
#'
#' Figure out gene locations to make room for gene names.
#' Written original by Dan Gatti 2013-02-13
#'
#' @param locs tbl of gene information
#' @param xlim X axis limits
#' @param text_size size of text (default 3)
#' @param str_rect character spacing on left and right of rectangles (default c("iW","i"))
#' @param n_rows desired number of rows (default 10)
#' @param plot_width width of default plot window (in inches)
#'
#' @return list object used by \code{\link{gene_plot}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#'    Daniel Gatti, \email{Dan.Gatti@@jax.org}
#' @references \url{https://github.com/dmgatti/DOQTL/blob/master/R/gene.plot.R}
#' @keywords utilities
get.gene.locations = function(locs, xlim, text_size=3,
                              str_rect=c("iW","i"),
                              n_rows=10, plot_width=6, ...) {
  
  ## Goal. Get rid of reassignment.
  ## Tighten code.
  ## Translate position to mm. Relate to gglot2 default text size of 7.
  ##
  
  ## This is in mm. Really want it in xlim coordinates.
  ## Should perhaps depend on text size.
  str_width <- function(chars, text_fudge = 4) {
    strwidth(chars, "inches") * diff(xlim) * text_size /
      (text_fudge * plot_width)
  }
  str_rect <- rep(str_rect, length.out=2)
  str_i <- str_width(str_rect[2])
  str_iW <- str_width(str_rect[1])
  str_name <- str_width(locs$Name)
  
  ## Line genes up sequentially in columns.
  
  nrows <- 1 # Number of rows in the current plot.
  row <- n_rows  # Number of rows that we would like.
  ymin <- 1  # Lowest y-value for the gene closest to the bottom of the plot.
  iter <- 0  # Number of iterations.
  
  # We need at least enough rows in the plot to fit the data.
  while(nrows < row & iter < 20) {
    nrows = n_rows
    row = 1
    x = min(locs$start)
    n_locs <- nrow(locs)
    
    ## Try to fill in the genes without collisions.
    
    ## Set up pointer to use instead of extra vector.
    ptr <- seq_len(n_locs)
    bottom <- rep(1L, n_locs)
    keep <- rep(TRUE, n_locs)
    
    i = 1    # Gene counter.
    while(i <= n_locs & any(keep)) {
      # Find gene with minimum start past x.
      wh = which(locs$start[keep] >= x)
      
      # If we found one, give it a bottom row and remove from keep.
      if(length(wh)) {
        # idx indexes the current gene row in tmp.
        idx = seq_len(n_locs)[keep][min(wh)]
        ptr[i] <- idx
        bottom[i] <- row
        keep[idx] <- FALSE
        
        if(any(keep)) {
          # Update x. If past plot, advance row and reset x.
          x <- locs$stop[idx] + str_name[idx] + str_iW
          if(x > xlim[2]) {
            row <- row + 1
            x <- min(locs$start[keep])
          }
        }
        i <- i + 1
      } else {
        ## If no gene past current X position, advance to next row
        ## and reset X position to left edge of plot.
        row <- row + 1
        x <- min(locs$start[keep])
      }
    }
    iter <- iter + 1
  }
  bottom[ptr] <- bottom
  list(nrows = nrows, row = row, bottom = bottom, nudge = str_i)
}
