#' Set up colors for patterns or points
#'
#' @param scan1output output of linear mixed model for \code{phename} (see \code{\link[qtl2]{scan1}})
#' @param snpinfo Data frame with snp information
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}.
#' @param col Color of other points, or colors for patterns
#' @param pattern allele pattern as of form \code{AB:CDEFGH}
#' @param show_all_snps show all SNPs if \code{TRUE}
#' @param col_hilit Color of highlighted points
#' @param drop_hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#' @param maxlod Maximum LOD for drop of \code{drop_hilit}
#'
#' @return list of \code{col} and \code{pattern}.
#'
#' @importFrom dplyr desc distinct filter group_by mutate summarize ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#'
color_patterns_set <- function(scan1output, snpinfo, patterns,
                               col, pattern, show_all_snps,
                               col_hilit, drop_hilit, maxlod) {
  if(!show_all_snps) { # reduce to sdp for distinct SNPs
    distinct_snps <- match(rownames(scan1output), snpinfo$snp)
  }
  if(patterns != "none") {
    # col rank-ordered by decreasing lod for pheno and pattern
    if(is.na(drop_hilit) || is.null(drop_hilit))
      drop_hilit <- 1.5
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      pattern <- pattern[distinct_snps]
    }
    # Group by phenotype and pattern to find groups within drop_hilit of maxlod.
    # Get at most 7 distinct hi groups in order of decreasing lodGpPhe.
    upattern <- dplyr::distinct(
      dplyr::filter(dplyr::arrange(
        dplyr::ungroup(
          dplyr::summarize(
            dplyr::group_by(
              # Pivot longer LODs by phenotype across patterns.
              tidyr::pivot_longer(
                dplyr::mutate(
                  as.data.frame(unclass(scan1output)),
                  pattern = pattern),
                -pattern, names_to = "pheno", values_to = "lod"),
              .data$pheno, .data$pattern),
            lodPhenoPattern = max(.data$lod),
            hi = (.data$lodPhenoPattern >= maxlod - drop_hilit))),
        dplyr::desc(.data$lodPhenoPattern)),
        .data$hi),
      pattern)$pattern
    
    # Allow at most 8 patterns to be retained.
    upattern <- subset(upattern, seq_along(upattern) < 8)
    lpattern <- sort(unique(pattern))

    if(missing(col) || length(col) != length(lpattern)) {
      col <- match(lpattern, upattern, nomatch = 8)
      names(col) <- ifelse(col == 8, "other", lpattern)
    }
  } else { # patterns == "none"
    # pattern is pheno-specific indication of below or above drop_hilit threshold
    # col set for pheno and pattern
    nphe <- dim(scan1output)[2]
    col <- rep(col, len = nphe)
    names(col) <- dimnames(scan1output)[[2]]
    # Highlight above drop_hilit?
    if(!is.na(drop_hilit) && !is.null(drop_hilit)) {
      pattern <- nphe * (scan1output >= maxlod - drop_hilit) + col(scan1output)
      col <- c(col, rep(col_hilit, len = nphe))
      names(col) <- seq_along(col)
    } else {
      pattern <- NULL
    }
  }

  # Shape parameter for point
  if(is.null(shape <- snpinfo$type)) {
    shape <- "snp"
  } else {
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      shape <- shape[distinct_snps]
    }
  }

  list(pattern = pattern, col = col, shape = shape)
}

#' Set up col, pattern, shape and group for plotting.
#'
#' @param scan1ggdata data frame to be used for plotting
#' @param lod matrix of LOD scores by position and pheno
#' @param pattern allele pattern of form \code{AB:CDEFGH}
#' @param col Color for \code{color} column in \code{scan1ggdata}
#' @param shape Shape for \code{shape} column in \code{scan1ggdata}
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}
#' @param facet use \code{\link[ggplot2]{facet_wrap}} if not \code{NULL}
#'
#' @return data frame \code{scan1ggdata} with additional objects.
#'
#' @importFrom dplyr filter mutate rename
color_patterns_pheno <- function(scan1ggdata,
                                 lod,
                                 pattern,
                                 col,
                                 shape,
                                 patterns,
                                 facet = NULL) {
  # Modify columns in scan1ggdata for plotting.
  #   col = pheno if patterns == "none"
  #         pattern otherwise
  #   group = chr_pheno for discrete line drawing
  #   facets = pattern if patterns == "none"
  #           pheno otherwise

  # If there is only one pheno, then pattern becomes pheno.
  if(!is.null(pattern)) {
    # If provided, pattern has to be same size as lod.
    labels <- NULL
    if(!is.matrix(lod))
      lod <- matrix(lod, length(lod), 1)
    if(!is.matrix(pattern)) {
      if(is.factor(pattern))
        labels <- levels(pattern)
      pattern <- matrix(pattern, length(pattern), 1)
    } else {
      if(length(pattern) == nrow(lod)) {
        pattern <- matrix(pattern, length(pattern), 1)
      }
    }
    if(!(nrow(pattern) == nrow(lod))) {
      stop("pattern must have same length as lod")
    }
    pattern_list <- color_patterns_other(pattern, lod, col, labels)
    if(ncol(lod) == 1) {
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   group = paste(.data$chr, pattern_list$pattern, sep = "_"))
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   color = pattern_list$color)
    } else {
      # *** Need to do following after collapsing colors.
      # *** But make sure the "pattern" remains intact for connecting lines.
      # Set up col and facet columns.
      if(is.null(facet) | facet %in% c("pheno","geno")) { # col=pattern, facet=pheno
        scan1ggdata <- dplyr::rename(scan1ggdata,
                                     facets = .data$pheno)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(.data$chr, pattern_list$pattern, sep = "_"))
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     color = pattern_list$color)
      } else { # facet == "pattern": col=pheno, facet=pattern
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     facets = pattern_list$color)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(.data$chr, .data$pheno, sep = "_"))
        scan1ggdata <- dplyr::rename(scan1ggdata,
                                     color = .data$pheno)
      }
    }
  } else { # no pattern provided
    if(!is.null(facet)) {
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   facets = .data$pheno)
    }
    ## group makes chr and col combination distinct for plotting.
    scan1ggdata <- dplyr::mutate(scan1ggdata,
                                 group = paste(.data$chr, .data$pheno, sep = "_"))
    ## If want facet, assume it is pheno.
    scan1ggdata <- dplyr::rename(scan1ggdata,
                                 color = .data$pheno)
  }
  # Make sure group is ordered.
  scan1ggdata <- dplyr::mutate(scan1ggdata,
                               group = ordered(.data$group, levels = unique(.data$group)))
  
  # shape for plotting
  if(is.null(shape)) {
    shape <- "snp"
  }

  scan1ggdata <- dplyr::mutate(scan1ggdata,
                               shape = rep(shape,
                                           length = nrow(scan1ggdata)))

  if(patterns == "hilit") {
    scan1ggdata <- dplyr::filter(scan1ggdata,
                                 pattern_list$color != pattern_list$other)
  }

  scan1ggdata
}

color_patterns_other <- function(pattern, lod, col,
                                 labels = NULL) {

  # Extend pattern if needed to have same length as lod.
  pattern <- c(matrix(pattern, nrow(lod), ncol(lod)))

  # Order levels of pattern by descinding lod
  lod <- rep(lod, len = length(pattern))
  if(is.null(labels))
    labels <- unique(pattern[order(-lod)])
  pattern <- factor(pattern, levels = labels)

  ## Reduce pattern to levels based on names(col) if provided
  other <- "other"
  color <- pattern
  if(!is.null(col)) {
    col <- rep(col, length = length(labels))
    # Names of colors may be subset of labels.
    if(!is.null(names(col))) {
      # The names(col) should be a subset of labels, with additional "other" name.
      # If some names of col agree with pheno, then collapse pheno.
      # Primarily used for patterns with plot_snpasso.
      m <- match(labels, names(col), nomatch = 0)
      if(!all(m == 0)) {
        if(any(m == 0)) {
          other <- names(col[-m])
          if(all(other == other[1])) {
            # If all other are the same, then change those color labels to other.
            tmpfn <- function(color, m, labels, other) {
              tmp <- as.character(color)
              tmp[tmp %in% labels[m == 0]] <- other[1]
              factor(tmp, levels = c(labels[m>0], other[1]))
            }
            color = tmpfn(pattern, m, labels, other[1])
          }
        }
      }
    }
  }
  list(pattern = pattern, color = color, other = other[1])
}

#' Set up col, pattern and group for plotting.
#'
#' @param scan1ggdata data frame to be used for plotting
#' @param col Color for \code{color} column in \code{scan1ggdata}
#' @param palette for colors (default \code{NULL} uses \code{"Dark2"} from \code{RColorBrewer} package)
#'
#' @return list of \code{colors} and \code{shapes}.
#'
#' @importFrom dplyr filter mutate rename
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
color_patterns_get <- function(scan1ggdata, col, palette=NULL) {
  # Set up colors using palette.
  labels <- levels(scan1ggdata$color)
  if(is.null(col)) {
    col <- seq_along(labels)
  }
  if(is.null(names(col))) {
    # Assume col without names match labels more or less.
    col <- rep(col, length.out = min(length(col), length(labels)))
    names(col) <- rep(labels, length.out = length(col))
  }
  if(!all(labels %in% names(col))) {
    col <- seq_along(labels)
    names(col) <- labels
  } else {
    col <- col[labels]
  }

  if(is.numeric(col)) {
    if(is.null(palette)) palette <- "Dark2"
    colors <-
      RColorBrewer::brewer.pal(
        RColorBrewer::brewer.pal.info[palette,"maxcolors"], palette)
    ncolors <- length(colors)
    col <- colors[1 + ((col-1) %% ncolors)]
  }

  shapes <- c(snp=1,indel=23,SV=25,INS=25,DEL=24,INV=22)
  if(is.null(shape <- scan1ggdata$shape)) {
    shapes <- shapes[1]
  } else {
    shape <- factor(shape)
    ## See http://sape.inf.usi.ch/quick-reference/ggplot2/shape
    ## Add diamond shape to any overlooked above.
    tmp <- levels(shape) %in% names(shapes)
    if(any(!tmp)) {
      newshapes <- levels(shape)[!tmp]
      shapes <- c(shapes, rep(21,length(newshapes)))
      names(shapes)[-(seq_along(newshapes))] <- newshapes
    }
    shapes <- shapes[levels(shape)]
  }

  list(colors = col, shapes = shapes)
}
