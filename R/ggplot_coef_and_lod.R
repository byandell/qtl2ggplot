# two-panel plot with both coefficients and LOD scores
# (for a single chromosome)
#
# calls ggplot_coef and ggplot_scan1
# internal function that is called by plot_coef
ggplot_coef_and_lod <-
    function(x, map, columns=NULL, col=NULL, scan1_output,
             gap=25, ylim=NULL, bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects",
             ylab_lod="LOD score", ylim_lod=NULL, col_lod="slateblue",
             vlines=NULL, main=FALSE,
             maxpos = NULL, maxcol = 1,
             legend.position = "right",
             top_panel_prop=0.65,
             lodcolumn = 1,
             facet = NULL, facet_lod = NULL,
             pattern = NULL, 
             pattern_lod = pattern,
             legend.position_lod = legend.position,
             ...)
{
    if(is.null(map)) stop("map is NULL")
      if(!is.list(map)) map <- list(" " = map) # if a vector, treat it as a list with no names
      
    # align scan1 output and map
    tmp <- qtl2::align_scan1_map(x, map)
    x <- tmp$scan1
    map <- tmp$map
    
    if(nrow(x) != length(unlist(map)))
      stop("nrow(x) [", nrow(x), "] != number of positions in map [",
           length(unlist(map)), "]")
    
    # also, match markers and use map in coefficients object
    # this seems clumsy and does not work well for multiple traits
    mar_in_coef <- rownames(x)
    mar_in_scan1 <- rownames(scan1_output)
    
    ## subset individuals in scan1 output.
    wh <- which(mar_in_scan1 %in% mar_in_coef)
    scan1_output <- modify_object(scan1_output, 
                                  scan1_output[wh, , drop=FALSE])

    ## also fix pattern
    if(!is.null(pattern)) {
      pattern <- pattern[wh,, drop=FALSE]
      rownames(pattern) <- rownames(scan1_output)
    }
    mis_mar <- !(mar_in_coef %in% mar_in_scan1)
    if(any(mis_mar)) {
        n_new <- sum(mis_mar)
        new_lod <- matrix(NA, nrow=sum(mis_mar), ncol=ncol(scan1_output))
        rownames(new_lod) <- mar_in_coef[mis_mar]
        scan1_output <- rbind(scan1_output, new_lod)[mar_in_coef,]

        ## fix pattern with missing values
        if(!is.null(pattern)) {
          pattern <- rbind(pattern, new_lod)[mar_in_coef,]
          rownames(pattern) <- NULL
        }
    }

    if(is.null(maxpos)) { # include vertical line at max lod
      maxpos <- qtl2::max_scan1(scan1_output, map, lodcolumn = 1)$pos[1]
    }

    # 2 x 1 panels
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow = 2,
                                   heights=c(top_panel_prop,
                                             1-top_panel_prop))))

    p1 <- ggplot_coef(x, map=map, columns=columns, col=col, scan1_output=NULL,
                      add=FALSE, gap=gap, ylim=ylim, bgcolor=bgcolor,
                      altbgcolor=altbgcolor, ylab=ylab,
                      vines=vlines, main=main,
                      legend.position = legend.position,
                      maxpos = maxpos, maxcol = maxcol,
                      facet = facet, pattern = pattern, ...) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())

    p2 <- ggplot_scan1(scan1_output, map=map, lodcolumn=lodcolumn, col=col_lod,
                     gap=gap, vines = vlines,
                     legend.position = legend.position_lod,
                     pattern = pattern_lod, facet = facet_lod, ...) +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank())

    if(!is.na(maxpos))
      p2 <- p2 + ggplot2::geom_vline(xintercept = maxpos,
                                     linetype=2,
                                     col = maxcol)
    print(p1, 
          vp = grid::viewport(layout.pos.row = 1,
                              layout.pos.col = 1))
    print(p2,
          vp = grid::viewport(layout.pos.row = 2,
                              layout.pos.col = 1))
}
