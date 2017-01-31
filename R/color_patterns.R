# Rethinking color and group
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
#   pheno group=1,2 maxval hi col

# could also make col character as color name or hex.

# set up colors for patterns or points
color_patterns_set <- function(scan1output, patterns,
                               col, group, show_all_snps,
                               col.hilit, drop.hilit, maxlod) {
  if(patterns != "none") {
    if(is.na(drop.hilit) || is.null(drop.hilit))
      drop.hilit <- 1.5
    if(!show_all_snps) { # reduce to sdp for distinct SNPs
      tmp <- match(scan1output$map[[1]],
                   scan1output$snpinfo[[1]]$pos_Mbp)
      group <- group[tmp]
    }
    # Set color for all SDPs with max below maxlod-drop.hilit to color 8.
    # Now modifying for multiple phenotypes. Gets tricky.
    lod <- as.data.frame(scan1output$lod)
    lod$group <- group
    lod <- tidyr::gather(lod, pheno, lod, -group)
    group_hi <- dplyr::arrange(
      dplyr::ungroup(
        dplyr::summarize(
          dplyr::group_by(lod, pheno, group),
          maxval = max(lod),
          hi = maxval >= maxlod - drop.hilit)),
      dplyr::desc(maxval))

    if(missing(col) || length(col) != nrow(group_hi)) {
      group_hi <- dplyr::ungroup(
        dplyr::mutate(
          dplyr::group_by(group_hi, pheno),
          col = pmax(ifelse(hi, rank(-maxval), 8), 8)))
      col <- group_hi$col
    }
    names(col) <- names(group_hi)
    names(col)[!group_hi] <- "other"
  } else {
    # Highlight above drop.hilit?
    if(!is.na(drop.hilit) && !is.null(drop.hilit)) {
      group <- (1:2)[(scan1output$lod >= maxlod-drop.hilit)+1]
    } else {
      group <- NULL
    }
    col <- c(col, col.hilit)
  }
  list(group = group, col = col)
}

color_patterns_pheno <- function(scan1ggdata, col, patterns) {
  labels <- levels(scan1ggdata$pheno)
  if(!is.null(col)) {
    col <- rep(col, length = length(labels))
    # Names of colors may be subset of labels.
    if(!is.null(names(col))) {
      # The names(col) should be a subset of labels, with additional "other" name.
      # If some names of col agree with pheno, then collapse pheno.
      # Primarily used for patterns with plot_snpasso.
      m <- match(names(col), labels, nomatch = 0)
      if(!all(m == 0)) {
        if(any(m == 0)) {
          other <- names(col[m == 0])
          if(all(other == other[1])) {
            tmp <- as.character(scan1ggdata$pheno)
            tmp[tmp %in% labels[m == 0]] <- other[1]
            scan1ggdata$pheno <- factor(tmp, c(labels[m], other[1]))
          }
        }
      }
    }
  }
  if(patterns == "hilit") {
    scan1ggdata <- dplyr::filter(scan1ggdata, pheno != other[1])
  }

  scan1ggdata
}
color_patterns_get <- function(scan1ggdata, col, palette=NULL) {
  # Set up colors using palette.
  labels <- levels(scan1ggdata$pheno)
  if(!is.null(col)) {
    if(!is.null(names(col)))
      col <- col[labels]
    names(col) <- NULL
  } else {
    col <- seq_along(labels)
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
