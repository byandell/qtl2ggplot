#' Plot a genome scan
#'
#' Plot LOD curves for a genome scan
#'
#' @param map Map of pseudomarker locations.
#' @param lod Matrix of lod (or other) values.
#' @param gap Gap between chromosomes.
#' @param col Colors for points or lines, with labels.
#' @param shape Shapes for points. 
#' @param pattern Use to group values for plotting (default = \code{NULL}); typically provided by \code{\link{plot_snpasso}} internal routine.
#' @param facet Plot facets if multiple phenotypes and group provided (default = \code{NULL}).
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}.
#'
#' @param ... Additional graphics parameters.
#'
#' @importFrom ggplot2 ggplot aes xlim ylim xlab ylab ggtitle
#' facet_grid facet_wrap geom_line geom_point theme geom_rect
#' scale_x_continuous coord_cartesian
#' theme element_rect element_blank
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate rename
#' @importFrom stringr str_replace
#' @importFrom rlang .data
#' @rdname ggplot_scan1
#'
ggplot_scan1_internal <-
  function(map, lod, gap = 25,
           col=NULL,
           shape=NULL,
           pattern = NULL, facet = NULL,
           patterns = c("none","all","hilit"),
           ...)
  {
    patterns <- match.arg(patterns)
    scan1ggdata <- make_scan1ggdata(map, lod, gap, col, pattern, shape,
                                    facet, patterns)

    ## Make sure we don't invoke facets if no column present.
    if(!match("facets", names(scan1ggdata), nomatch = 0))
      facet <- NULL
    
    ggplot_scan1_create(map, gap, col, shape, scan1ggdata, facet, ...)
  }

make_scan1ggdata <- function(map, lod, gap, col, pattern, shape,
                             facet, patterns) {
  # set up chr and xpos with gap.
  xpos <- unlist(map) # map_to_xpos(map, gap)
  chr <- factor(rep(names(map), sapply(map, length)),
                levels = names(map))

  # make data frame for ggplot
  rownames(lod) <- NULL # make sure duplicates do not mess us up for multiple traits
  # Make sure colnames of lod are unique for pivot_longer. 
  tmp <- colnames(lod)
  colnames(lod) <- paste0(letters[seq_along(tmp)], tmp)
  scan1ggdata <- data.frame(xpos=xpos, chr=chr, lod,
                            check.names = FALSE)
  scan1ggdata <- tidyr::pivot_longer(scan1ggdata, 
                               -(1:2), names_to = "pheno", values_to = "lod")
  scan1ggdata <- dplyr::mutate(scan1ggdata, 
                               pheno = stringr::str_replace(.data$pheno, "^[a-z]", ""))
  # make sure order of pheno is preserved.
  scan1ggdata <- dplyr::mutate(scan1ggdata,
                               pheno = ordered(.data$pheno, levels = unique(.data$pheno)))

  ## facet if more than one pheno or set by user.
  if(ncol(lod) > 1 & !is.null(pattern)) {
    # If facet is not NULL, pattern  column of scan1ggdata is used to facet.
    # That column is either pheno or pattern, set in color_patterns_pheno.
    if(is.null(facet))
      facet <- "pheno"
  }
  ## Set up col, group and (optional) facet in scan1ggdata.
  ## Column pheno becomes either col or facet
  color_patterns_pheno(scan1ggdata,
                       lod, pattern, col, shape,
                       patterns, facet)
}

ggplot_scan1_create <-
  function(map, gap, col, shape, scan1ggdata, facet,
           bgcolor, altbgcolor,
           lwd=1,
           pch = col_shape$shapes,
           cex=1,
           point_fill = "transparent",
           xlab=NULL, ylab="LOD score",
           xaxt = ifelse(onechr, "y", "n"), 
           yaxt = "y",
           palette = "Dark2",
           xlim=NULL, ylim=NULL, main=FALSE,
           legend.position =
             ifelse(length(levels(scan1ggdata$color)) == 1, "none", "right"),
           legend.title="pheno",
           lines=TRUE, points=!lines,
           scales = c("free_x","free"),
           ...)
  {

    scales <- match.arg(scales)
    
    # Extra arguments
    onechr <- (length(map)==1) # single chromosome

    if(is.null(xlab)) {
      xlab <- "Position"
    }

    if(onechr & !is.null(xlim)) {
      scan1ggdata <- dplyr::filter(scan1ggdata,
                                   .data$xpos >= xlim[1] & .data$xpos <= xlim[2])
      if(!nrow(scan1ggdata)) {
        warning(paste("no plot data in range", xlim[1], "to", xlim[2]))
        return(NULL)
      }
    }
    
    # make ggplot aesthetic with limits and labels
    p <- ggplot2::ggplot(scan1ggdata,
                         ggplot2::aes(x = .data$xpos, y = .data$lod,
                                      col = .data$color,
                                      shape = shape,
                                      group = .data$group)) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
    
    # gap between chromosomes
    p <- p +
      ggplot2::theme(panel.spacing = grid::unit(gap / 10000, "npc"))

    # Facets (if multiple phenotypes and groups).
    if(all(levels(scan1ggdata$chr) == " ")) {
      if(!is.null(facet)) {
        p <- p + ggplot2::facet_wrap( ~ facets, scales = scales)
      }
    } else {
      if(!is.null(facet)) {
        p <- p + ggplot2::facet_grid(facets ~ chr, scales = scales, space = "free")
      } else {
        p <- p + ggplot2::facet_grid( ~ chr, scales = scales, space = "free")
      }
    }

    # color palette, point shapes and legend titles
    col_shape <- color_patterns_get(scan1ggdata, col, palette)
    p <- p +
      ggplot2::scale_color_manual(name = legend.title,
                                  values = col_shape$colors)
    p <- p +
      ggplot2::scale_shape_manual(name = "SV Type",
                                  labels = names(pch),
                                  values = pch)

    # add legend if requested
    p <- p +
      ggplot2::theme(legend.position = legend.position)

    # include axis labels?
    if(yaxt == "n") {
      p <- p + 
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank())
    }
    # X axis
    if(xaxt == "n") {
      p <- p + 
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(), 
          axis.ticks.x = ggplot2::element_blank())
    }
    if(!onechr)
      xlim <- NULL
    p <- p + 
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

    # remove y axis?
    if(yaxt == "n") {
      p <- p + 
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank())
    }
    # grid lines
    p <- ggplot_grid_lines(p, onechr, ...)
    # add box just in case
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black",
                                             fill=NA))

    # add main as title if provided
    # or use name from lod if only one column
    if(!is.logical(main)) {
      title <- main
      main <- TRUE
    }
    if(main) {
      if(title == "") {
        # create title from groups if only 1
        group_names <- levels(scan1ggdata$group)
        if(length(group_names) == 1) {
          p <- p +
            ggplot2::ggtitle(group_names) +
            ggplot2::theme(legend.position = "none")
        }
      } else {
        p <- p +
          ggplot2::ggtitle(title)
      }
    }

    ## Add lines and/or points.
    if(lines) {
      p <- p + ggplot2::geom_line(size = lwd)
    }
    if(points) {
      p <- p + ggplot2::geom_point(size = cex,
                                   fill = point_fill)
    }

    p
  }
