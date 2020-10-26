### R/qtl2ggplot

[Karl Broman](https://kbroman.org) & [Brian Yandell](http://pages.stat.wisc.edu/~yandell/)

R/qtl2ggplot is a reimplementation of [qtl2plot](https://github.com/rqtl/qtl2plot) for data visualization. It includes all functions in qtl2plot, but now using [ggplot2](https://ggplot2.tidyverse.org/) and related routines. See
[R/qtl2](https://kbroman.org/qtl2/) (aka qtl2) for the bigger story of the qtl2 suite of routines.

---

### Installation

First install [R/qtl2](https://kbroman.org/qtl2/) and a set of package
dependencies, including the
[devtools](https://github.com/r-lib/devtools) package.

    install.packages(c("tidyverse", "RColorBrewer", "grid", "qtl2"))

Then, install qtl2ggplot using `devtools::install_github()`.

    library(devtools)
    install_github("byandell/qtl2ggplot")

To install vignettes:

    install_github("byandell/qtl2ggplot", build_vignettes = TRUE)

---

### Vignettes

- [qtl2ggplot](https://github.com/byandell/qtl2ggplot/blob/master/vignettes/qtl2ggplot.Rmd)

---

#### License

[Licensed](License.md) under [GPL-3](https://www.r-project.org/Licenses/GPL-3).
