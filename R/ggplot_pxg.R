#' Plot phenotype vs genotype
#'
#' Plot phenotype vs genotype for a single putative QTL and a single phenotype.
#'
#' @param geno Vector of genotypes, as produced by
#' \code{\link[qtl2]{maxmarg}} with specific \code{chr} and
#' \code{pos}.
#' @param pheno Vector of phenotypes.
#' @param sort If TRUE, sort genotypes from largest to smallest.
#' @param SEmult If specified, interval estimates of the within-group
#' averages will be displayed, as \code{mean +/- SE * SEmult}.
#' @param pooledSD If TRUE and \code{SEmult} is specified, calculated
#' a pooled within-group SD. Otherwise, get separate estimates of
#' the within-group SD for each group.
#' @param jitter Amount to jitter the points horizontally, if a vector
#' of length > 0, it is taken to be the actual jitter amounts
#' (with values between -0.5 and 0.5).
#' @param seg_width Width of segments at the estimated within-group averages
#' @param seg_lwd Line width used to plot estimated within-group averages
#' @param seg_col Line color used to plot estimated within-group averages
#' @param bgcolor Background color for the plot.
#' @param hlines Locations of horizontal grid lines.
#' @param hlines_col Color of horizontal grid lines
#' @param hlines_lty Line type of horizontal grid lines
#' @param hlines_lwd Line width of horizontal grid lines
#' @param vlines_col Color of vertical grid lines
#' @param vlines_lty Line type of vertical grid lines
#' @param vlines_lwd Line width of vertical grid lines
#' @param force_labels If TRUE, force all genotype labels to be shown.
#' @param alternate_labels If TRUE, place genotype labels in two rows
#' @param omit_points If TRUE, omit the points, just plotting the averages (and, potentially, the +/- SE intervals).
#' @param ... Additional graphics parameters, passed to \code{\link[graphics]{plot}}.
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#'
#' @export
#' @importFrom ggplot2 aes element_rect geom_jitter geom_point geom_segment geom_vline ggplot theme
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom stats lm runif sd
#' @importFrom rlang .data
#'
#' @seealso \code{\link{plot_coef}}
#'
#' @examples
#' # load qtl2 package for data and genoprob calculation
#' library(qtl2)
#'
#' # read data
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#'
#' # inferred genotype at a 28.6 cM on chr 16
#' geno <- maxmarg(probs, map, chr=16, pos=28.6, return_char=TRUE)
#'
#' # plot phenotype vs genotype
#' ggplot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)))
#'
#' # include +/- 2 SE intervals
#' ggplot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          SEmult=2)
#'
#' # plot just the means
#' ggplot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          omit_points=TRUE)
#'
#' # plot just the means +/- 2 SEs
#' ggplot_pxg(geno, log10(iron$pheno[,1]), ylab=expression(log[10](Liver)),
#'          omit_points=TRUE, SEmult=2)
ggplot_pxg <-
  function(geno, pheno, sort=TRUE, SEmult=NULL, pooledSD=TRUE,
           jitter=0.2, bgcolor="gray90",
           seg_width=0.4, seg_lwd=2, seg_col="black",
           hlines=NULL, hlines_col="white", hlines_lty=1, hlines_lwd=1,
           vlines_col="gray80", vlines_lty=1, vlines_lwd=3,
           force_labels=TRUE, alternate_labels=FALSE,
           omit_points=FALSE, ...)
  {
    if(length(jitter) > 1) {
      stop("length(jitter) > 1 not allowed")
    }
    if(length(jitter) == 1) {
      if(jitter < 0 || jitter > 0.5)
        stop("jitter should be in [0,0.5]")
    }

    ggplot_pxg_internal <-
      function(geno, pheno, jitter = 0.2, bgcolor="gray90",
               seg_width=0.2, seg_lwd=2, seg_col="slateblue",
               hlines=NULL, hlines_col="white", hlines_lty=1, hlines_lwd=1,
               vlines_col="gray80", vlines_lty=1, vlines_lwd=3,
               xlim=c(0.5, length(unique(geno))+0.5), ylim=NULL,
               xaxs="i", yaxs="r", xlab="Genotype", ylab="Phenotype",
               mgp=c(2.6, 0.3, 0), mgp.x=mgp, mgp.y=mgp, las=1,
               pch=21, bg="blue",
               alpha=0.25, ...)
      {
        ids <- qtl2::get_common_ids(geno, pheno)
        dat <- data.frame(
          geno = geno[ids],
          pheno = pheno[ids])
        
        datu <- mean_pxg(geno, pheno, dat)
        datu$index <- seq_len(nrow(datu))
        
        p <- ggplot2::ggplot(dat) +
          ggplot2::aes(x=geno, y=pheno) +
          ggplot2::geom_point( # placeholder to set up x axis and data frame
            ggplot2::aes(x=geno, y=mean),
            data = datu,
            col = "transparent",
            inherit.aes = FALSE) +
          ggplot2::geom_vline(
            xintercept = datu$index,
            col = vlines_col, 
            size = vlines_lwd / 2)

        if(!omit_points) {
          p <- p +
            ggplot2::geom_jitter(width = jitter, height = 0, alpha = alpha,
                                 col = bg)
        }

        if(seg_width > 0) {
          seg_width <- seg_width / 2
          p <- p +           
            ggplot2::geom_segment(
              ggplot2::aes(x    = .data$index - seg_width,
                           xend = .data$index + seg_width,
                           y    = .data$mean,
                           yend = .data$mean),
              data = datu,
              size = seg_lwd / 2,
              col = seg_col,
              inherit.aes = FALSE)
          
          if(!is.null(SEmult)) {
            p <- p + 
              ggplot2::geom_segment(
                ggplot2::aes(x    = .data$index,
                             xend = .data$index,
                             y    = .data$mean - .data$se * SEmult,
                             yend = .data$mean + .data$se * SEmult),
                data = datu,
                size = seg_lwd / 2,
                col = seg_col,
                inherit.aes = FALSE)
          }
        }
          
        # Add box for each chr.
        p <- p +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(colour = "black",
                                                 fill=NA))
        
        p
      }
    
    ggplot_pxg_internal(geno, pheno, jitter, bgcolor=bgcolor,
                      seg_width=seg_width, seg_lwd=seg_lwd, seg_col=seg_col,
                      hlines=hlines, hlines_col=hlines_col, hlines_lty=hlines_lty, hlines_lwd=hlines_lwd,
                      vlines_col=vlines_col, vlines_lty=vlines_lty, vlines_lwd=vlines_lwd, ...)
  }
#' Mean phenotype by genotype
#' 
#' @param geno Vector of genotypes, as produced by
#' \code{\link[qtl2]{maxmarg}} with specific \code{chr} and
#' \code{pos}.
#' @param pheno Vector of phenotypes.
#' @param dataframe Supplied data frame, or constructed from \code{geno} and \code{pheno} if \code{NULL}.
#' 
#' @rdname ggplot_pxg
#' @export
#' 
mean_pxg <- function(geno, pheno, dataframe = NULL) {
  if(is.null(dataframe)) {
    ids <- qtl2::get_common_ids(geno, pheno)
    dataframe <- data.frame(
      geno = geno[ids],
      pheno = pheno[ids])
  }
  
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(dataframe, geno),
      mean = mean(pheno, na.rm = TRUE),
      sd = sd(pheno, na.rm = TRUE),
      se = sd / sqrt(sum(!is.na(pheno)))))
}