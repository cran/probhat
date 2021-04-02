\name{70_moment-based_statistics}
\alias{ph.mean}
\alias{ph.sd}
\alias{ph.var}
\alias{ph.skewness}
\alias{ph.kurtosis}
\alias{moment}
\alias{central.moment}
\alias{standardized.moment}
\alias{raw.moment}
\title{Moment-Based Statistics}
\description{Compute moment-based statistics from probability distributions.}
\usage{
ph.mean (sf, \dots, n.intervals=200)
ph.sd (sf, \dots, n.intervals=200)
ph.var (sf, \dots, n.intervals=200)
ph.skewness (sf, \dots, n.intervals=200)
ph.kurtosis (sf, \dots, n.intervals=200)

moment (sf, nth, \dots, n.intervals=200)

central.moment (sf, nth, \dots, n.intervals=200)
standardized.moment (sf, nth, \dots, n.intervals=200)

raw.moment (sf, nth, about=0, \dots, n.intervals=200)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this is a univariate PMF or spline-based CDF.\cr\cr Refer to the references and see also sections.}
\item{nth}{Integer, the nth moment}.
\item{about}{Numeric, the about constant for raw moments.}
\item{n.intervals}{Integer.\cr In the discrete case, ignored.\cr In the continuous case, the number of intervals, used in the numerical approximation.}
\item{\dots}{Ignored.}
}
\details{
The mean/sd/var/skewness/kurtosis functions all call the moment function.

If the moment function is called with nth equal zero, it returns one.\cr
If called with nth=1 (the mean), it computes the the first raw moment.\cr
If called with nth=2 (the variance), it computes the second central moment.\cr
If called with nth=3 (the skewness), it computes the third standardized moment.\cr
If called with nth=4 (the kurtosis), it computes the fourth standardized moment.\cr
And if called with nth>4, it also computes the nth standardized moment.

Note that currently, the standard deviation, variance and higher moments, should should not be regarded as accurate.
}
\value{A single numeric value.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Empirical-Like Distributions}

\link{probmv}, \link{rng}

\link{quartiles}, \link{ntiles}, \link{ph.median}, \link{ph.quantile}\cr
\link{ph.mode}, \link{ph.modes}
}
\examples{
prep.ph.data ()

cFh <- cdfuv.cks (height)

ph.mean (cFh)
ph.sd (cFh)
ph.skewness (cFh)
ph.kurtosis (cFh)
}
