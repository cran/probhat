\name{15_multivariate_probabilities}
\alias{probmv}
\title{Multivariate Probabilities}
\description{Compute multivariate probabilities, from multivariate CDFs.}
\usage{
probmv (F, a, b)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{F}{A cdfmv.cks object.}
\item{a}{A vector (or matrix) of lower limits, corresponding to each variable. If a is a matrix, then each row relates to one region and each column relates one variable.}
\item{b}{Same as a, except upper limits.}
}
\examples{
data.prep ()

#multivariate cdf
Fh = cdfmv.cks (trees)

#approximate first and third quartiles
a = c (28, 22, 0.55)
b = c (38, 24, 1.05)

#multivariate probability
probmv (Fh, a, b)
}
