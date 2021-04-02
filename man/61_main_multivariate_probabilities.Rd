\name{61_main_multivariate_probabilities}
\alias{pwith}
\alias{pwith.cksuv}
\alias{pwith.cksmv}
\title{Multivariate Probabilities}
\description{
Compute the probability of observing a point within a rectangular/cuboidal/etc region.

Currently, these functions are limited to nonconditional \link{CKS} models.\cr
The function \link{probmv} can be used for almost-arbitrary continuous CDFs, but is deprecated. 
}
\usage{
pwith (\dots)

\method{pwith}{cksuv}(sf, xlim = c (a, b), \dots, a=-Inf, b=Inf)
\method{pwith}{cksmv}(sf, xlim = cbind (a, b), \dots, a=-Inf, b=Inf)
}
\arguments{
\item{sf}{
	A suitable function object.\cr
	(Here, this is a cksuv or cksmv object, which are nonconditional \link{CKS} objects.)\cr\cr
	PDFs, CDFs and quantile functions, are all allowed.}
\item{xlim}{
	A length-two numeric vector or a two-column numeric matrix.\cr
	A standard vector is the same as a one-row matrix.\cr
	And the number of rows needs to match the number of random variables.\cr\cr
	This represents the lower and upper limits, corresponding to the limits of integration, if one was integrating the corresponding PDF.\cr\cr
	By default, this is just the combination a and b, meaning that limits can be specified in two ways.}
\item{a, b}{Numeric vectors, giving the lower and upper limits, respectively.}
\item{\dots}{Ignored.}
}
\details{
Refer to the vignette for more information.

Note that evaluating a spline-based CDF is likely to faster, ignoring the construction of the object.
}
\value{
A single numeric value.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{pdist}, \link{probmv}

\link{cdfmv.cks}, \link{cdfmvc.cks}

\link{ph.mean}\cr
\link{moment}

\link{rng}

\link{quartiles}, \link{ntiles}, \link{ph.median}, \link{ph.quantile}\cr
\link{ph.mode}, \link{ph.modes}\
}
\examples{
prep.ph.data ()

cfh3 <- pdfmv.cks (trees2)

xlim <- matrix (c (
    22, 24,    #height in 22 to 24
    28, 38,    #girth  in 28 to 38
    0.55, 1.05 #volume in 0.55 to 1.05
    ),, 2, byrow=TRUE, dimnames = list (colnames (trees), c ("a", "b") ) )

pwith (cfh3, xlim=xlim)
}
