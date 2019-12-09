\name{17_parameter_estimates}
\alias{ph.mean}
\alias{ph.sd}
\alias{ph.var}
\alias{ph.skewness}
\alias{ph.kurtosis}
\alias{ph.median}
\alias{ph.quantile}
\alias{ph.mode}
\alias{ph.moment}
\alias{ph.moment.2}
\title{Parameter Estimates}
\description{Compute parameter (and parameter-like) estimates from univariate or conditional distributions.}
\usage{
#from PMF of continuous CDF
ph.mean (f, n=200)
ph.sd (f, n=200)
ph.var (f, n=200)
ph.skewness (f, n=200)
ph.kurtosis (f, n=200)
ph.moment (f, nth, n=200)
ph.moment.2 (f, nth, about, n=200)

#from quantile function
ph.median (F.inv)
ph.quantile (F.inv, p)

#from PMF or PDF
ph.mode (f)
}
\details{Refer to the vignette for more information.

The ph.mean/sd/var/skewness/kurtosis functions all call the ph.moment function.

The ph.median/quantile functions are the same as evaluating the quantile function.

The ph.mode function returns a single global mode. If there's more than one global mode, then the first is returned.

The ph.moment function (but not the ph.moment.2 function) automatically centers and sometimes standardizes the moment, depending on what moment is being computed. The zeroth moment returns one. The first moment (which gives the mean) is zero-centered and unstandardized. The second moment (which gives the variance) is mean-centered and unstandardized. Higher moments are mean-centered and standardized).

Note that currently, the standard deviation, variance and higher moments, are not reliable, and should only be used as as an exploratory tool, for the purpose of comparing different probability distributions.
}
\arguments{
\item{f}{For the mean, sd, variance and moments, either a PMF or continuous CDF. For the mode, either a PMF or PDF.}
\item{F.inv}{A quantile function.}
\item{p}{The probability.}
\item{nth}{The nth moment}.
\item{about}{.}
\item{n}{In the continuous case (otherwise ignored), the number of intervals, used in the numerical approximation.}
}
\examples{
data.prep ()

conditions = c (Girth=30, Height=22)
Fh = cdfc.cks (trees, conditions=conditions)
ph.mean (Fh)
}
