\name{44_other_plot_methods}
\alias{ph.plotf.ph4.gset}
\alias{ph.plotf.ph4.mset}
\title{Other Plot Methods}
\description{Plot methods for distribution sets.\cr\cr NOTE THAT THESE FUNCTIONS WILL BE REPLACED.}
\usage{
\method{ph.plotf}{ph4.gset}(sfs, \dots, span.win=FALSE, legend=TRUE, colors)
\method{ph.plotf}{ph4.mset}(sfs, \dots, nr, nc, colors)
}
\arguments{
\item{sfs}{A ph4.gset or ph4.mset object.}
\item{span.win}{Logical, if true, plots (most notably CDFs) extend across the entire plotting window.}
\item{legend}{Logical, if true, add a legend to the plot.}
\item{colors}{Character vector, colors for each distribution.}
\item{nr, nc}{Optional integers, number of panels.}
\item{\dots}{Ignored.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Distribution Sets}
}
\examples{
prep.ph.data ()

plot (ph4.pdfuv.gset.cks (species, sepal.length) )
plot (ph4.qfuv.mset.el (trees2), nr=2, nc=2)
}
