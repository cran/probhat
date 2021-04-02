\name{31_distribution_sets}
\alias{Distribution Sets}
\alias{cset}
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
\alias{pdfmv.gset.cks}
\title{Distribution Sets}
\description{
	Sets of distributions.\cr
	Note that distributions sets can also be produced by conditional \link{CKS} models.\cr\cr
	NOTE THAT THESE OBJECTS ARE SUBJECT TO CHANGE, AND HAVE HAD LIMITED TESTING.\cr
	NOTE THAT as.list SHOULD BE CALLED ON THE RESULTING OBJECTS, UNLESS USED WITH PLOTTING FUNCTIONS.}
\usage{
#categorical sets
#(grouped by a categorical variable)
pdfuv.gset.cks (g, x, \dots)
cdfuv.gset.cks (g, x, \dots)
qfuv.gset.cks (g, x, \dots)
cdfuv.gset.el (g, x, \dots)
qfuv.gset.el (g, x, \dots)

#marginal sets
pdfuv.mset.cks (x, \dots, bw, smoothness=1)
cdfuv.mset.cks (x, \dots, bw, smoothness=1)
qfuv.mset.cks (x, \dots, bw, smoothness=1)
cdfuv.mset.el (x, \dots)
qfuv.mset.el (x, \dots)

pdfmv.gset.cks (g, x, \dots)
}
\arguments{
\item{g}{A character vector, with the same length as x.\cr Or an object which can be coerced to such a vector.\cr If a list, its first element is used.}
\item{x}{Vector (for categorical sets) or matrix (for marginal sets).}
\item{bw, smoothness}{Bandwidth and smoothness parameters, same as pdfuv.cks, except that they can be an m-length vector, where m is the number of variables, equal to number of columns in x.}
\item{\dots}{Other arguments for the corresponding constructor.}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

These functions construct distribution set objects, which are lists of probability distributions.
}
\value{A ph4.gset or ph4.mset object.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{ph.plotf.ph4.gset}, \link{ph.plotf.ph4.mset}
}
\examples{
prep.ph.data ()

gs <- pdfuv.gset.cks (species, sepal.length)
ms <- qfuv.mset.el (trees2)

plot (gs)
plot (ms, nr=2, nc=2)

#distance matrix
pdist (gs)

#please wrap inside as.list, to use each object 
gs <- as.list (gs)
plot (gs [[1]])
}
