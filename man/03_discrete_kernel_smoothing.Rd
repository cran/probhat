\name{03_discrete_kernel_smoothing}
\alias{pmfuv.dks}
\alias{cdfuv.dks}
\alias{qfuv.dks}
\title{Models via Discrete Kernel Smoothing}
\description{
Fit (estimated) probability distributions, with discrete kernel smoothing, from data.

Note that the default bandwidth (of one) results in an (unsmoothed) frequency distribution. Also note that only unique x values are retained, and frequencies are automatically aggregated, if there are any duplicates.
}
\usage{
pmfuv.dks (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
cdfuv.dks (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
qfuv.dks (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
}
\details{
Refer to the vignette for more information.

Note that x, h or both is needed.

If x is supplied and h is omitted, h defaults to the frequencies of x.
However, if h is supplied and x is omitted, then x defaults to 1:n.
}
\value{
Self-referencing functions of the form:

function (x) = ...\cr
function (y) = ...

Where x is an integer-valued vector, and y is a numeric vector, with values in the interval [0, 1].
}
\arguments{
\item{x}{Integer-valued vector or a single-column integer-valued matrix, preferably with a column name, and optionally with row names.}
\item{h}{Vector of counts, which can be fractional.}
\item{sequence}{Ignored.}
\item{kernel}{Constructor for a discrete kernel object, or the kernel object itself.}
\item{bw}{A positive odd integer.}
\item{lower}{Optional integer, the lower truncation value. If missing, defaults to the smallest x value. Setting lower to -Inf (negative infinity), prevents truncation.}
}
\examples{
data.prep ()

fh = pmfuv.dks (traffic.x, traffic.h, bw=25, lower=0)
fh (20)
}
