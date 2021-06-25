#' List of scan1coef objects
#'
#' Create a list of scan1coef objects using \code{\link[qtl2]{scan1coef}}.
#'
#' @param phe data frame of phenotypes
#' @param probs genotype probabilities object for one chromosome from \code{\link[qtl2]{calc_genoprob}}
#' @param K list of length 1 with kinship matrix
#' @param covar matrix of covariates
#' @param blups Create BLUPs if \code{TRUE}
#' @param center Make intercept constant if \code{TRUE}
#' @param ... Additional arguments
#'
#' @return object of class \code{listof_scan1coef}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' # read data
#' iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#'
#' # insert pseudomarkers into map
#' map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
#'
#' # Ensure that covariates have names attribute
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#'
#' # Calculate scan1coef on all phenotypes,
#' # returning a list of \code{\link{scan1coef}} objects
#' out <- listof_scan1coef(probs[,7], iron$pheno, addcovar = covar, center = TRUE)
#' 
#' # Plot coefficients for all phenotypes
#' ggplot2::autoplot(out, map[7], columns = 1:3)
#' 
#' # Summary of coefficients at scan peak
#' scan_pr <- qtl2::scan1(probs[,7], iron$pheno)
#' summary(out, scan_pr, map[7])
#' 
#' @export
#' @importFrom qtl2 scan1coef scan1blup
#' 
listof_scan1coef <- function(probs, phe, K=NULL, covar=NULL, blups = FALSE,
                             center = FALSE, ...) {
  eff <- list()
  phename <- dimnames(phe)[[2]]
  scan1fn <- ifelse(blups, 
                  qtl2::scan1blup,
                  qtl2::scan1coef)

  for(pheno in phename) {
    out <- scan1fn(probs, phe[, pheno, drop=FALSE], K, covar, ...)
    if(center) {
      intcol <- match("intercept", colnames(out))
      mall <- mean(out[, intcol], na.rm = TRUE)
      nall <- ncol(out) - 1
      out[, -intcol] <- apply(out[,-intcol], 2, function(x, y) x - mall + y, out[,intcol])
      out[, intcol] <- mall
    }
    eff[[pheno]] <- out
  }
  class(eff) <- c("listof_scan1coef", class(eff))
  eff
}

#' Summary of object of class listof_scan1coef
#'
#' Summary of object of class \code{\link{listof_scan1coef}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coef}
#' @param scan1_object object from \code{scan1}
#'
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2]{insert_pseudomarkers}}.
#'
#' @param coef_names names of effect coefficients (default is all coefficient names)
#' @param center center coefficients if \code{TRUE}
#' @param ... arguments for \code{\link[qtl2]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
#' @importFrom dplyr bind_cols
summary_listof_scan1coef <-
  function(object, scan1_object, map,
           coef_names = dimnames(object[[1]])[[2]],
           center = TRUE,
           ...) {
  phename <- names(object)
  chr_id <- names(map)

  # Summary of phenos on chr; reorder and rename by phename
  sum_chr <- summary(scan1_object, map, lodcolumn=phename, chr=chr_id)
  m <- match(phename, sum_chr$pheno)
  sum_chr <- sum_chr[m,]
  sum_chr$pheno <- phename
  
  pos <- sum_chr$pos
  names(pos) <- phename
  wh <- apply(scan1_object,2,function(x) which.max(x)[1])
  sum_coef <- as.data.frame(matrix(NA,
                                   nrow(sum_chr),
                                   length(coef_names),
                                   dimnames=list(NULL,coef_names)))
  for(pheno_id in phename) {
    if(!is.null(object[[pheno_id]])) {
      ## Need to save these in data frame.
      tmp <- object[[pheno_id]][wh[pheno_id], seq_along(coef_names)]
      if(center)
        tmp <- tmp - mean(tmp)
      sum_coef[match(pheno_id, sum_chr$pheno),] <- tmp
    }
  }
  dplyr::bind_cols(sum_chr,sum_coef)
}
#' @method summary listof_scan1coef
#' @rdname listof_scan1coef
#' @export
summary.listof_scan1coef <- function(object, ...)
  summary_listof_scan1coef(object, ...)

#' Summary of object of class listof_scan1coef
#'
#' Summary of object of class \code{\link{listof_scan1coef}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coef}
#' @param scan1_object object from \code{scan1}
#' @param ... arguments for \code{\link[qtl2]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
summary_scan1coef <-
  function(object, scan1_object, map, ...) {
    if(!inherits(object, "listof_scan1coef")) {
      object <- list(object)
      names(object) <- dimnames(scan1_object)[[2]][1]
    }
    summary_listof_scan1coef(object, scan1_object, map, ...)
  }

#' @method summary scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export summary.scan1coef
#'
summary.scan1coef <- function(object, ...)
  summary_scan1coef(object, ...)

#' Subset of object of class listof_scan1coef
#'
#' Subset of object of class \code{\link{listof_scan1coef}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param x object of class \code{listof_scan1coef}
#' @param elements indexes or names of list elements in x
#' @param ... ignored
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
subset_listof_scan1coef <- function(x, elements, ...) {
  x_class <- class(x)
  class(x) <- "list"
  x <- x[elements]
  class(x) <- x_class
  x
}

#' @method subset listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export subset.listof_scan1coef
#'
subset.listof_scan1coef <- function(x, ...)
  subset_listof_scan1coef(x, ...)

#' @method [ listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export [.listof_scan1coef
#'
`[.listof_scan1coef` <- function(x, ...)
  subset_listof_scan1coef(x, ...)
