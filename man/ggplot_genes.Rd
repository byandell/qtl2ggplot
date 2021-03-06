% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot_genes.R
\name{ggplot_genes}
\alias{ggplot_genes}
\alias{autoplot.genes}
\title{Plot gene locations for a genomic interval}
\usage{
ggplot_genes(
  genes,
  xlim = NULL,
  minrow = 4,
  padding = 0.2,
  colors = c("black", "red3", "green4", "blue3", "orange"),
  ...
)

\method{autoplot}{genes}(x, ...)
}
\arguments{
\item{genes}{Data frame containing \code{start} and \code{stop} in
bp, \code{strand} (as \code{"-"}, \code{"+"}, or \code{NA}), and
\code{Name}.}

\item{xlim}{x-axis limits (in Mbp)}

\item{minrow}{Minimum number of rows of genes}

\item{padding}{Proportion to pad with white space around the genes}

\item{colors}{Vectors of colors, used sequentially and then re-used.}

\item{...}{Optional arguments passed to \code{\link[graphics]{plot}}.}

\item{x}{Object of class \code{genes}}
}
\value{
None.
}
\description{
Plot gene locations for a genomic interval, as rectangles with gene
symbol (and arrow indicating strand/direction) below.
}
\examples{
filename <- file.path("https://raw.githubusercontent.com/rqtl",
                      "qtl2data/master/DOex", 
                      "c2_genes.rds")
tmpfile <- tempfile()
download.file(filename, tmpfile, quiet=TRUE)
gene_tbl <- readRDS(tmpfile)
unlink(tmpfile)

ggplot_genes(gene_tbl)

}
\keyword{hgraphics}
