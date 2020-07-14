\name{71_multivariate_probabilities}
\alias{probmv}
\title{Multivariate Probabilities}
\description{Compute probabilities, from multivariate CDFs.}
\usage{
probmv (sf, a, b)
}
\arguments{
\item{sf}{A cdfmv.cks or cdfmvc.cks object.}
\item{a, b}{Numeric vectors (or matrices) of lower and upper limits, corresponding to each variable.\cr If they're matrices, then each row defines the limits for one region and each column defines the limits for one variable.}
}
\details{Refer to the vignette for more information.}
\value{
A single numeric value (if a and b are both standard vectors) and a numeric vector (if either a or b are matrices).
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{cdfmv.cks}, \link{cdfmvc.cks}

\link{ph.mean}\cr
\link{moment}

\link{rng}

\link{quartiles}, \link{ntiles}, \link{ph.median}, \link{ph.quantile}\cr
\link{ph.mode}, \link{ph.modes}\
}
\examples{
ph.data.prep ()

cFh3 <- cdfmv.cks (trees)

q <- matrix (c (
    22, 24,    #height in 22 to 24
    28, 38,    #girth  in 28 to 38
    0.55, 1.05 #volume in 0.55 to 1.05
    ),, 2, byrow=TRUE, dimnames = list (colnames (trees), c ("a", "b") ) )

probmv (cFh3, q [,1], q [,2])
}
