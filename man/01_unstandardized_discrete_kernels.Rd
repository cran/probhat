\name{01_unstandardized_discrete_kernels}
\alias{binomial.kernel}
\title{Unstandardized Discrete Kernels}
\description{Unstandardized discrete kernels for discrete kernel smoothing.}
\usage{
binomial.kernel ()
}
\details{
This is a constructor for a binomial kernel object.

Currently, this object is a named list, however, it's possible that it may be changed to an S4 object, in the future.

This object has $pmf and $cdf components, each of which are functions, giving its PMF and CDF, respectively.

Here, the PMF is symmetric about zero (in contrast to how a binomial distribution is normally defined), and has zero-probability outside the interval [-hbw, hbw], where:\cr
hbw = (bw - 1) / 2
}
\value{
A named list with two component functions.

Each component function has the form:

function (bw, x) = ...

Where bw is a positive odd integer, and x is an integer vector.
}
\examples{
dkernel = binomial.kernel ()
dkernel$pmf (13, 0)
dkernel$cdf (13, 0)
}
