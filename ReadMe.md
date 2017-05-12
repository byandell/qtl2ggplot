### R/qtl2ggplot

[![Build Status](https://travis-ci.org/rqtl/qtl2plot.svg?branch=master)](https://travis-ci.org/rqtl/qtl2plot)

[Karl Broman](http://kbroman.org) & [Brian Yandell](http://www.stat.wisc.edu/~yandell)

R/qtl2ggplot is a reimplementation of [qtl2plot](https://github.com/rqtl/qtl2plot) (for data visualization). It includes all functions in qtl2plot, but now using [ggplot2](http://ggplot2.org/) and related routines. See
[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) for the bigger story of the qtl2 suite of routines.

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You first need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed)

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))
    
You will also need the following packages for qtl2ggplot:

    install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "grid"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github(paste0("rqtl/qtl2", c("geno", "scan")))

Once you have installed these, install qtl2ggplot as

    install_github(paste0("byandell/", c("CCSanger", "DOread")))
    install_github(paste0("byandell/qtl2", c("pattern", "ggplot")))

---

### Vignettes

coming ...

---

#### License

[Licensed](License.md) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).
