# two-panel plot with both coefficients and LOD scores
# (for a single chromosome)
#
# calls plot_coef and plot_scan1
# internal function that is called by plot_coef
#' @importFrom grid grid.layout grid.newpage pushViewport viewport
plot_coef_and_lod <-
    function(x, columns=NULL, col=NULL, scan1_output,
             gap=25, ylim=NULL, bgcolor="gray90", altbgcolor="gray85",
             ylab="QTL effects",
             ylab_lod="LOD score", ylim_lod=NULL, col_lod="slateblue",
             xaxt=NULL,
             vlines=NULL, main=FALSE,
             maxlod=TRUE, maxcol = 1,
             legend.position = "none",
             top_panel_prop=0.65, ...)
{
    # also, match markers and use map in coefficients object
    mar_in_coef <- rownames(x$coef)
    mar_in_scan1 <- rownames(scan1_output$lod)
    scan1_output$lod <- scan1_output$lod[mar_in_scan1 %in% mar_in_coef, , drop=FALSE]
    scan1_output$map <- list("1"=x$map)
    mis_mar <- !(mar_in_coef %in% mar_in_scan1)
    if(any(mis_mar)) {
        n_new <- sum(mis_mar)
        new_lod <- matrix(NA, nrow=sum(mis_mar), ncol=scan1_output$lod)
        rownames(new_lod) <- mar_in_coef[mis_mar]
        scan1_output$lod <- rbind(scan1_output$lod, new_lod)[mar_in_coef,]
    }
    
    if(maxlod) { # include vertical line at max lod
      maxpos <- max(scan1_output, lodcolumn = 1)$pos[1]
    }
    
    # 2 x 1 panels
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow = 2,
                                   heights=c(top_panel_prop, 
                                             1-top_panel_prop))))
    
    p1 <- plot_coef(x, columns=columns, col=col, scan1_output=NULL,
                    add=FALSE, gap=gap, ylim=ylim, bgcolor=bgcolor,
                    altbgcolor=altbgcolor, ylab=ylab,
                    xaxt="n", vines=vlines, main=main, legend.position = legend.position, ...)
    if(maxlod)
      p1 <- p1 + ggplot2::geom_vline(xintercept = maxpos, 
                                     linetype=2,
                                     col = maxcol)
    print(p1, vp = grid::viewport(layout.pos.row = 1,
                                  layout.pos.col = 1))

    p2 <- plot_scan1(scan1_output, lodcolumn=1, col=col_lod,
                     gap=gap, vines = vlines, legend.position = legend.position, ...)
    if(maxlod)
      p2 <- p2 + ggplot2::geom_vline(xintercept = maxpos, 
                                     linetype=2,
                                     col = maxcol)
    print(p2, vp = grid::viewport(layout.pos.row = 2,
                                  layout.pos.col = 1))
}
