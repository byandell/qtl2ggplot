# Rethinking color and pattern
# The group identifies all the possible groups
# Each group has a col, but multiple groups may have the same col.
# Groups and cols differ by phenotype, but have the same quality.
# Need to think if col is by rank order of group or by group identifier.
# (Think patterns.)
# So group could be the same size as lod (see ggplot_scan1 handling)
# but col could be a list the same length as number of phenos.

# The modified lod df below has pheno-group-lod connection.
# Do we want to keep group_hi, which has group-col assignment?
# Don't alter scan1 object, to keep compatible with Karl's approach.

# Need to take care of legacy (col, col.hilit).
# could create small group_hi set
#   pheno group=1,2 lodGpPhe hi col

# could also make col character as color name or hex.

# set up colors for patterns or points

#' @importFrom dplyr desc distinct filter group_by mutate summarize ungroup
#' @importFrom tidyr gather
color_patterns_set <- function(scan1output, lodcolumns, patterns,
                               col, pattern, show_all_snps,
                               col.hilit, drop.hilit, maxlod) {
  if(patterns != "none") {
    if(is.na(drop.hilit) || is.null(drop.hilit))
      drop.hilit <- 1.5
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      tmp <- match(scan1output$map[[1]],
                   scan1output$snpinfo[[1]]$pos_Mbp)
      pattern <- pattern[tmp]
    }
    # Set color for all SDPs with max below maxlod-drop.hilit to color 8.
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
            lodGpPhe = max(lod),
            hi = lodGpPhe >= maxlod - drop.hilit)),
        dplyr::desc(lodGpPhe)),
        hi),
      pattern)$pattern
    upattern <- subset(upattern, seq_along(upattern) < 8)
    lpattern <- sort(unique(pattern))
    
    if(missing(col) || length(col) != length(lpattern)) {
      col <- match(lpattern, upattern, nomatch = 8)
      names(col) <- ifelse(col == 8, "other", lpattern)
    }
  } else { # patterns == "none"
    nphe <- dim(scan1output$lod)[2]
    col <- rep(col, len = nphe)
    # Highlight above drop.hilit?
    if(!is.na(drop.hilit) && !is.null(drop.hilit)) {
      pattern <- nphe * (scan1output$lod >= maxlod - drop.hilit) + col(scan1output$lod)
      col <- c(col, rep(col.hilit, len = nphe))
    } else {
      pattern <- NULL
    }
  }
  list(pattern = pattern, col = col)
}

#' Set up pheno and pattern, and chr_pheno, for colors.
#' 
#' @importFrom tidyr gather
#' @importFrom dplyr filter mutate rename
color_patterns_pheno <- function(scan1ggdata, 
                                 lod, 
                                 pattern, 
                                 col, 
                                 patterns, 
                                 facet_var) {
  # Modify columns in scan1ggdata for plotting.
  #   col = pheno if patterns == "none"
  #         pattern otherwise
  #   group = chr_pheno for discrete line drawing
  #   facet = pattern if patterns == "none"
  #           pheno otherwise

  # If there is only one pheno, then pattern becomes pheno.
  if(!is.null(pattern)) {
    # If provided, pattern has to be same size as lod.
    if(!is.matrix(pattern))
      pattern <- matrix(pattern, length(pattern), 1)
    if(!is.matrix(lod))
      lod <- matrix(lod, length(lod), 1)
    if(!(nrow(pattern) == nrow(lod))) {
      stop("pattern must have same length as lod")
    }
    pattern_list <- color_patterns_other(pattern, lod, col)
    if(ncol(lod) == 1) {
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   group = paste(chr, pattern_list$pattern, sep = "_"))
      scan1ggdata <- dplyr::mutate(scan1ggdata,
                                   color = pattern_list$color)
    } else {
      # *** Need to do following after collapsing colors.
      # *** But make sure the "pattern" remains intact for connecting lines.
      # Set up col and facet columns.
      facet_set <- facet_var
      if(is.null(facet_var) | facet_var == "pheno") { # col=pattern, facet=pheno
        scan1ggdata <- dplyr::rename(scan1ggdata, 
                                     facet = pheno)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(chr, pattern_list$pattern, sep = "_"))
        scan1ggdata <- dplyr::mutate(scan1ggdata, 
                                     color = pattern_list$color)
      } else { # facet_var == "pattern": col=pheno, facet=pattern
        scan1ggdata <- dplyr::mutate(scan1ggdata, 
                                     facet = pattern_list$color)
        scan1ggdata <- dplyr::mutate(scan1ggdata,
                                     group = paste(chr, pheno, sep = "_"))
        scan1ggdata <- dplyr::rename(scan1ggdata, 
                                     color = pheno)
      }
    }
  } else { # no pattern provided
    if(!is.null(facet_var)) {
      scan1ggdata <- dplyr::mutate(scan1ggdata, 
                                   facet = pheno)
    }
    ## group makes chr and col combination distinct for plotting.
    scan1ggdata <- dplyr::mutate(scan1ggdata,
                                 group = paste(chr, pheno, sep = "_"))
    ## If want facet, assume it is pheno.
    scan1ggdata <- dplyr::rename(scan1ggdata, 
                                 color = pheno)
  }
  
  if(patterns == "hilit") {
    scan1ggdata <- dplyr::filter(scan1ggdata, 
                                 pattern_list$color != pattern_list$other)
  }
  
  scan1ggdata
}

color_patterns_other <- function(pattern, lod, col) {
  
  # Extend pattern if needed to have same length as lod.
  pattern_df <- as.data.frame(matrix(pattern, nrow(lod), ncol(lod)))
  names(pattern_df) <- dimnames(lod)[[2]]
  pattern_df <- tidyr::gather(pattern_df, pheno, pattern)
  pattern <- pattern_df$pattern
  
  # Order levels of pattern by descinding lod
  lod <- rep(lod, len = length(pattern))
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

color_patterns_get <- function(scan1ggdata, col, palette=NULL) {
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
  col
}
