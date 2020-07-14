\name{27_distribution_sets}
\alias{Distribution Sets}
\alias{pdfuv.gset.cks}
\alias{cdfuv.gset.cks}
\alias{qfuv.gset.cks}
\alias{cdfuv.gset.el}
\alias{qfuv.gset.el}
\alias{pdfuv.mset.cks}
\alias{cdfuv.mset.cks}
\alias{qfuv.mset.cks}
\alias{cdfuv.mset.el}
\alias{qfuv.mset.el}
\title{Distribution Sets}
\description{Sets of distributions.\cr\cr THESE FUNCTIONS HAVE HAD LIMITED TESTING.}
\usage{
#categorical sets
#(grouped by a categorical variable)
pdfuv.gset.cks (x, \dots, group.by)
cdfuv.gset.cks (x, \dots, group.by)
qfuv.gset.cks (x, \dots, group.by)
cdfuv.gset.el (x, \dots, group.by)
qfuv.gset.el (x, \dots, group.by)

#marginal sets
pdfuv.mset.cks (x, \dots, bw, smoothness=1)
cdfuv.mset.cks (x, \dots, bw, smoothness=1)
qfuv.mset.cks (x, \dots, bw, smoothness=1)
cdfuv.mset.el (x, \dots)
qfuv.mset.el (x, \dots)
}
\arguments{
\item{x}{Vector (for categorical sets) or matrix (for marginal sets).}
\item{group.by}{A character vector, with the same length as x.\cr Or an object which can be coerced to such a vector.\cr If a list, its first element is used.}
\item{bw, smoothness}{Bandwidth and smoothness parameters, same as pdfuv.cks, except that they can be an m-length vector, where m is the number of variables, equal to number of columns in x.}
\item{\dots}{Other arguments for the corresponding univariate constructor.}
}
\details{
THESE FUNCTIONS HAVE HAD LIMITED TESTED.

These functions construct distribution set objects, which are lists of probability distributions.
}
\value{A ph3.gset or ph3.mset object.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{plot.ph3.gset}, \link{plot.ph3.mset}
}
\examples{
ph.data.prep ()

plot (pdfuv.gset.cks (sepal.length, group.by=species) )
plot (qfuv.mset.el (trees), nr=2, nc=2)
}
