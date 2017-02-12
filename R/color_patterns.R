#' Set up colors for patterns or points
#'
#' @param scan1_output output of linear mixed model for \code{phename} (see \code{\link[qtl2scan]{scan1}})
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}.
#' @param col Color of other points, or colors for patterns
#' @param pattern allele pattern determined by \code{\link[CCSanger]{sdp_to_pattern}}
#' @param col.hilit Color of highlighted points
#' @param drop.hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#' @param show_all_snps show all SNPs if \code{TRUE}
#'
#' @return list of \code{col} and \code{pattern}
#'
#' @export
#' @importFrom dplyr desc distinct filter group_by mutate summarize ungroup
#' @importFrom tidyr gather
#'
color_patterns_set <- function(scan1output, patterns,
                               col, pattern, show_all_snps,
                               col.hilit, drop.hilit, maxlod) {
  if(!show_all_snps) { # reduce to sdp for distinct SNPs
    distinct_snps <- match(scan1output$map[[1]],
                 scan1output$snpinfo[[1]]$pos_Mbp)
  }
  if(patterns != "none") {
    # col rank-ordered by decreasing lod for pheno and pattern
    if(is.na(drop.hilit) || is.null(drop.hilit))
      drop.hilit <- 1.5
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      pattern <- pattern[distinct_snps]
    }
    # Group by phenotype and pattern to find groups within drop.hilit of maxlod.
    # Get at most 7 distinct hi groups in order of decreasing lodGpPhe.
    upattern <- dplyr::distinct(
      dplyr::filter(dplyr::arrange(
        dplyr::ungroup(
          dplyr::summarize(
            dplyr::group_by(
              # Gather LODs by phenotype across patterns.
              tidyr::gather(
                dplyr::mutate(
                  as.data.frame(scan1output$lod),
                  pattern = pattern),
                pheno, lod, -pattern),
              pheno, pattern),
            lodPhenoPattern = max(lod),
            hi = (lodPhenoPattern >= maxlod - drop.hilit))),
        dplyr::desc(lodPhenoPattern)),
        hi),
      pattern)$pattern
    upattern <- subset(upattern, seq_along(upattern) < 8)
    lpattern <- sort(unique(pattern))

    if(missing(col) || length(col) != length(lpattern)) {
      col <- match(lpattern, upattern, nomatch = 8)
      names(col) <- ifelse(col == 8, "other", lpattern)
    }
  } else { # patterns == "none"
    # pattern is pheno-specific indication of below or above drop.hilit threshold
    # col set for pheno and pattern
    nphe <- dim(scan1output$lod)[2]
    col <- rep(col, len = nphe)
    names(col) <- dimnames(scan1output$lod)[[2]]
    # Highlight above drop.hilit?
    if(!is.na(drop.hilit) && !is.null(drop.hilit)) {
      pattern <- nphe * (scan1output$lod >= maxlod - drop.hilit) + col(scan1output$lod)
      col <- c(col, rep(col.hilit, len = nphe))
      names(col) <- seq_along(col)
    } else {
      pattern <- NULL
    }
  }

  # Shape parameter for point
  if(is.null(type <- scan1output$snpinfo[[1]]$type)) {
    shape <- "SNP"
  } else {
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      type <- type[distinct_snps]
    }
    shape <- stringr::str_sub(type, 1, 3)
    shape[shape %in% c("InD","Ind")] <- "indel"
  }

  list(pattern = pattern, col = col, shape = shape)
}

#' Set up col, pattern, shape and group for plotting.
#'
#' @param scan1ggdata data frame to be used for plotting
#' @param lod matrix of LOD scores by position and pheno
#' @param pattern allele pattern determined by \code{\link[CCSanger]{sdp_to_pattern}}
#' @param col Color for \code{color} column in \code{scan1ggdata}
#' @param patterns Connect SDP patterns: one of \code{c("none","all","hilit")}
#' @param facet use \code{\link[ggplot2]{facet_wrap}} if not \code{NULL}
#'
#' @export
#' @importFrom tidyr gather
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
    if(!is.matrix(pattern)) {
      if(is.factor(pattern))
        labels <- levels(pattern)
      pattern <- matrix(pattern, length(pattern), 1)
    }
    if(!is.matrix(lod))
      lod <- matrix(lod, length(lod), 1)
    if(!(nrow(pattern) == nrow(lod))) {
      stop("pattern must have same length as lod")
    }
    pattern_list <- color_patterns_other(pattern, lod, col, labels)
    if(ncol(lod) == 1) {
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   group = paste(chr, pattern_list$pattern, sep = "_"))
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   color = pattern_list$color)
    } else {
      # *** Need to do following after collapsing colors.
      # *** But make sure the "pattern" remains intact for connecting lines.
      # Set up col and facet columns.
      if(is.null(facet) | facet %in% c("pheno","geno")) { # col=pattern, facet=pheno
        scan1ggdata <- dplyr::rename(scan1ggdata,
                                     facets = pheno)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(chr, pattern_list$pattern, sep = "_"))
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     color = pattern_list$color)
      } else { # facet == "pattern": col=pheno, facet=pattern
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     facets = pattern_list$color)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(chr, pheno, sep = "_"))
        scan1ggdata <- dplyr::rename(scan1ggdata,
                                     color = pheno)
      }
    }
  } else { # no pattern provided
    if(!is.null(facet)) {
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   facets = pheno)
    }
    ## group makes chr and col combination distinct for plotting.
    scan1ggdata <- dplyr::mutate(scan1ggdata,
                                 group = paste(chr, pheno, sep = "_"))
    ## If want facet, assume it is pheno.
    scan1ggdata <- dplyr::rename(scan1ggdata,
                                 color = pheno)
  }

  # shape for plotting
  if(is.null(shape)) {
    shape <- "SNP"
  }

  scan1ggdata <- dplyr::mutate(scan1ggdata,
                               shape = shape)

  if(patterns == "hilit") {
    scan1ggdata <- dplyr::filter(scan1ggdata,
                                 pattern_list$color != pattern_list$other)
  }

  scan1ggdata
}

color_patterns_other <- function(pattern, lod, col,
                                 labels = NULL) {

  # Extend pattern if needed to have same length as lod.
  pattern_df <- as.data.frame(matrix(pattern, nrow(lod), ncol(lod)))
  names(pattern_df) <- dimnames(lod)[[2]]
  pattern_df <- tidyr::gather(pattern_df, pheno, pattern)
  pattern <- pattern_df$pattern

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
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr filter mutate rename
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
color_patterns_get <- function(scan1ggdata, col, palette=NULL, shape) {
  # Set up colors using palette.
  labels <- levels(scan1ggdata$color)
  if(is.null(col)) {
    col <- seq_along(labels)
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

  shape <- factor(shape)
  ## See http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  ## Add diamond shape to any overlooked above.
  shapes <- c(SNP=96,indel=23,INS=25,DEL=24,INV=22)
  tmp <- levels(shape) %in% names(shapes)
  if(any(!tmp)) {
    newshapes <- levels(shape)[!tmp]
    shapes <- c(shapes, rep(21,length(newshapes)))
    names(shapes)[-(1:5)] <- newshapes
  }
  shapes <- shapes[levels(shape)]

  list(colors = col, shapes = shapes)
}
