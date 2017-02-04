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
#   pheno group=1,2 lodGpPhe hi col

# could also make col character as color name or hex.

# set up colors for patterns or points

#' @importFrom dplyr desc distinct filter group_by mutate summarize ungroup
#' @importFrom tidyr gather
color_patterns_set <- function(scan1output, lodcolumns, patterns,
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
    # Group by phenotype and group to find groups within drop.hilit of maxlod.
    # Get at most 7 distinct hi groups in order of decreasing lodGpPhe.
    ugroup <- dplyr::distinct(
      dplyr::filter(dplyr::arrange(
        dplyr::ungroup(
          dplyr::summarize(
            dplyr::group_by(
              # Gather LODs by phenotype across groups.
              tidyr::gather(
                dplyr::mutate(
                  as.data.frame(scan1output$lod),
                  group = group),
                pheno, lod, -group), 
              pheno, group),
            lodGpPhe = max(lod),
            hi = lodGpPhe >= maxlod - drop.hilit)),
        dplyr::desc(lodGpPhe)),
        hi),
      group)$group
    ugroup <- subset(ugroup, seq_along(ugroup) < 8)
    
    lgroup <- levels(factor(group))
    if(missing(col) || length(col) != length(lgroup)) {
      col <- match(lgroup, ugroup, nomatch = 8)
      names(col) <- ifelse(col == 8, "other", lgroup)
    }
  } else { # patterns == "none"
    nphe <- dim(scan1output$lod)[2]
    col <- rep(col, len = nphe)
    # Highlight above drop.hilit?
    if(!is.na(drop.hilit) && !is.null(drop.hilit)) {
      group <- nphe * (scan1output$lod >= maxlod - drop.hilit) + col(scan1output$lod)
      col <- c(col, rep(col.hilit, len = nphe))
    } else {
      group <- NULL
    }
  }
  list(group = group, col = col)
}

#' Set up pheno and group, and chr_pheno, for colors.
#' 
#' @importFrom tidyr gather
#' @importFrom dplyr filter mutate rename
color_patterns_pheno <- function(scan1ggdata, 
                                 lod, 
                                 group, 
                                 col, 
                                 patterns, 
                                 facet_var) {
  # Modify columns in scan1ggdata for plotting.
  #   col = pheno if patterns == "none"
  #         group (= pattern from plot_snpasso) otherwise
  #   group = chr_pheno for discrete line drawing
  #   facet = group if patterns == "none"
  #           pheno otherwise

  # If there is only one pheno, then group becomes pheno.
  if(!is.null(group)) {
    # If provided, group has to be same size as lod.
    if(!(length(group) == nrow(lod))) {
      if(!(length(group) == length(lod)))
        stop("group must have same length as lod")
    }
    if(ncol(lod) == 1) {
      scan1ggdata$col <- factor(group)
    } else {
      # Extend group if needed to have same length as lod.
      group_df <- as.data.frame(matrix(group, nrow(lod), ncol(lod)))
      names(group_df) <- dimnames(lod)[[2]]
      group_df <- tidyr::gather(group_df, pheno, group)
      group <- factor(group_df$group)
      
      # Set up col and facet columns.
      if(facet_var == "pheno") {
        scan1ggdata <- dplyr::rename(scan1ggdata, facet = col)
        scan1ggdata <- dplyr::mutate(scan1ggdata, col = group)
      } else {
        scan1ggdata <- dplyr::mutate(scan1ggdata, facet = group)
      }
    }
  } else {
    ## If want facet, assume it is col (= pheno).
    if(!is.null(facet_var)) {
      scan1ggdata <- dplyr::mutate(scan1ggdata, facet = col)
    }
  }
  
  ## group makes chr and col combination distinct for plotting.
  scan1ggdata <- dplyr::mutate(scan1ggdata,
                               group = paste(chr, col, sep = "_"))
  
  ## Reduce pheno to levels based on names(col) if provided
  labels <- levels(scan1ggdata$col)
  other <- "other"
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
            tmpfn <- function(col, m, labels, other) {
              tmp <- as.character(col)
              tmp[tmp %in% labels[m == 0]] <- other[1]
              factor(tmp, c(labels[m], other))
            }
            scan1ggdata <- dplyr::mutate(scan1ggdata, 
                                         col = tmpfn(col, m, labels, other[1]))
          }
        }
      }
    }
  }
  if(patterns == "hilit") {
    scan1ggdata <- dplyr::filter(scan1ggdata, col != other[1])
  }

  scan1ggdata
}

color_patterns_get <- function(scan1ggdata, col, palette=NULL) {
  # Set up colors using palette.
  labels <- levels(scan1ggdata$col)
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
