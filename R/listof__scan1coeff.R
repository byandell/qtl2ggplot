#' List of scan1coef objects
#'
#' Create a list of scan1coef objects using \code{\link[qtl2scan]{scan1coef}}.
#'
#' @param phe data frame of phenotypes
#' @param probs genotype probabilities object for one chromosome from \code{\link[qtl2geno]{calc_genoprob}}
#' @param K list of length 1 with kinship matrix
#' @param covar matrix of covariates
#'
#' @return object of class \code{listof_scan1coeff}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{listof_scan1coef(probs, phe)}
#'
#' @export
#' @importFrom qtl2scan scan1coef
listof_scan1coef <- function(probs, phe, K=NULL, covar=NULL) {
  eff <- list()
  phename <- dimnames(phe)[[2]]
  for(pheno in phename)
    eff[[pheno]] <- qtl2scan::scan1coef(probs, phe[, pheno, drop=FALSE], K, covar)
  class(eff) <- c("listof_scan1coef", class(eff))
  eff
}

#' Plot of object of class listof_scan1coeff
#'
#' Plot object of class \code{listof_scan1coeff}, which is a list of objects of class \code{scan1coef}.
#'
#' @param x object of class \code{listof_scan1coeff}
#' @param facet Plot facets if multiple phenotypes and group provided (default = \code{"pattern"}).
#' @param ... arguments for \code{\link[qtl2scan]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @method plot listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @importFrom dplyr bind_rows
plot.listof_scan1coef <- function(x, columns=NULL, col=NULL,
                                    facet = "pattern",
                                    ...) {
  # stretch out map
  map <- unlist(lapply(x, function(x) x$map))

  # bind together coef matrices
  coefs <- dplyr::bind_rows(lapply(x, function(x) as.data.frame(x$coef)),
                            .id = "pheno")
  pattern <- coefs$pheno
  coefs <- as.matrix(coefs[,-1])
  rownames(coefs) <- names(map)

  # Reform as one scan1coef object
  x <- x[[1]]
  x$map <- map
  x$coef <- coefs
  
  plot(x, columns, col, pattern = pattern,
       patterns = "all", facet = facet, ...)
}
#' Summary of object of class listof_scan1coeff
#'
#' Summary of object of class \code{\link{listof_scan1coeff}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coeff}
#' @param scan1_object object from \code{scan1}
#' @param coef_names names of effect coefficients (default is all coefficient names)
#' @param ... arguments for \code{\link[qtl2scan]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @importFrom dplyr bind_cols
summary.listof_scan1coef <-
  function(object, scan1_object,
           coef_names = dimnames(object[[1]]$coef)[[2]],
           ...) {
  phename <- names(object)
  chr_id <- names(scan1_object$map)
  sum_chr <- summary(scan1_object, phename, chr_id)
  pos <- sum_chr$pos[match(phename, sum_chr$pheno)]
  names(pos) <- phename
  wh <- apply(scan1_object$lod,2,function(x) which.max(x)[1])
  sum_coef <- as.data.frame(matrix(NA,
                                   nrow(sum_chr),
                                   length(coef_names),
                                   dimnames=list(NULL,coef_names)))
  for(pheno_id in phename) {
    if(!is.null(object[[pheno_id]])) {
      ## Need to save these in data frame.
      sum_coef[match(pheno_id, sum_chr$pheno),] <-
        object[[pheno_id]]$coef[wh[pheno_id],seq_along(coef_names)]
    }
  }
  dplyr::bind_cols(sum_chr,sum_coef)
}
#' Summary of object of class listof_scan1coeff
#'
#' Summary of object of class \code{\link{listof_scan1coeff}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coeff}
#' @param scan1_object object from \code{scan1}
#' @param ... arguments for \code{\link[qtl2scan]{plot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary scan1coef
#' @rdname listof_scan1coef
#' @export
summary.scan1coef <-
  function(object, scan1_object, ...) {
    object <- list(object)
    names(object) <- dimnames(scan1_object$lod)[[2]][1]
    summary.listof_scan1coef(object, scan1_object, ...)
  }
           