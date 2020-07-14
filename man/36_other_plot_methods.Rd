\name{36_other_plot_methods}
\alias{plot.ph3.gset}
\alias{plot.ph3.mset}
\title{Other Plot Methods}
\description{Plot methods for distribution sets.\cr\cr NOTE THAT THESE FUNCTIONS WILL BE REPLACED.}
\usage{
\method{plot}{ph3.gset}(x, \dots, legend=TRUE, colors)
\method{plot}{ph3.mset}(x, \dots, nr, nc, colors)
}
\arguments{
\item{x}{A ph3.gset or ph3.mset object.}
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
ph.data.prep ()

plot (pdfuv.gset.cks (sepal.length, group.by=species) )
plot (qfuv.mset.el (trees), nr=2, nc=2)
}
