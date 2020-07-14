\name{11_discrete_kernels}
\alias{Discrete Kernels}
\alias{binomial.dkernel}
\title{Discrete Kernels}
\description{Discrete kernels (without standardized intervals), for discrete kernel smoothing.

NOTE THAT THESE MAY BE REMOVED.}
\usage{
binomial.dkernel (bw)
}
\arguments{
\item{bw}{Integer, the bandwidth.\cr This should be a positive odd number.\cr (Even numbers are rounded up).}
}
\details{
NOTE THAT THESE MAY BE REMOVED.

This is a constructor for a binomial kernel object.

Currently, this object is a named list (with pmf and cdf components).

Here, the PMF is symmetric about zero (in contrast to how a binomial distribution is normally defined), and has zero-probability outside the interval [-hbw, hbw], where:\cr
hbw = (bw - 1) / 2

This object is used inside discrete kernel smoothing.\cr
(In which case, you provide the constructor, not the object).
}
\value{A dkernel object.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Continuous Kernels}

\link{plot.dkernel}

\link{Discrete Kernel Smoothing}
}
\examples{
k <- binomial.dkernel (5)
plot (k)

k$pmf (-2:2)
}
