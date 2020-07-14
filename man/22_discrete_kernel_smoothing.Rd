\name{22_discrete_kernel_smoothing}
\alias{Discrete Kernel Smoothing}
\alias{pmfuv.dks}
\alias{cdfuv.dks}
\alias{qfuv.dks}
\title{Discrete Kernel Smoothing Models}
\description{Fit probability distributions, via discrete kernel smoothing over integer-indexed frequency data.}
\usage{
pmfuv.dks (x, h, \dots,
    bw.method="ph.default",
    kernel=binomial.dkernel,
    bounded = c (TRUE, FALSE), freq=FALSE,
    lower, upper,
    bw, smoothness=1)

cdfuv.dks (x, h, \dots,
    bw.method="ph.default",
    kernel=binomial.dkernel,
    bounded = c (TRUE, FALSE),
    lower, upper,
    bw, smoothness=1)

qfuv.dks (x, h, \dots,
    bw.method="ph.default",
    kernel=binomial.dkernel,
    bounded = c (TRUE, FALSE),
    lower, upper,
    bw, smoothness=1)
}
\arguments{
\item{x}{Integer vector of integer-indexed bins.\cr
	Also, can be a single-column integer matrix, preferably with a column (variable) name, and optionally with row (bin) names.\cr
	In principle, this is required.\cr
	(It can be omitted if h supplied, but generates a warning).}
\item{h}{Optional numeric vector of frequencies (or weights), which can be fractional.\cr Defaults to a vector of ones.}
\item{bw.method}{String, the bandwidth selection method.\cr Refer to \link{Bandwidth Selection}.}
\item{kernel}{Constructor for a dkernel (discrete kernel) object.}
\item{bounded}{Logical vector of length one or two, determining whether the probability distribution is lower or upper bounded.}
\item{freq}{Logical, if true, the resulting function object returns frequencies, by default.}
\item{bw, smoothness}{Integer and numeric values, the bandwidth and smoothness parameters, respectively.\cr
	The bandwidth parameter should be a positive odd number.\cr (Even numbers are incremented).\cr
	If the bandwidth is missing, it's computed using bw.method (see above) and the smoothness.}
\item{lower, upper}{Optional integers, the lower and upper truncation values.\cr
	They default to the lowest and highest bin values.\cr
	Note that x (the bin values) must not be outside the truncation values.\cr
	If so, and error is produced.\cr\cr
	Ignored, unless the corresponding bounded values are true.}
\item{\dots}{Ignored.}
}
\details{
Refer to the vignette for more information.

Note that if x has non-unique values, then duplicated x (and their h) values are aggregated.\cr
And currently, any row names will be ignored.
}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernels}

\link{Succinct Constructors}\cr
\link{Continuous Kernel Smoothing}, \link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{is.dks}, \link{print.phmodel}, \link{plot.dksuv}

\link{Bandwidth Selection}
}
\examples{
ph.data.prep ()

dFht <- qfuv.dks (traffic.bins, traffic.freq, lower=0)
dFht (0.5)
}
