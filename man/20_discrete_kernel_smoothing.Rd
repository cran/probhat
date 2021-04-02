\name{20_discrete_kernel_smoothing}
\alias{DKS}
\alias{Discrete Kernel Smoothing}
\alias{pmfuv.dks}
\alias{cdfuv.dks}
\alias{qfuv.dks}
\title{Discrete Kernel Smoothing Models}
\description{
Fit probability distributions, via discrete kernel smoothing over integer-indexed frequency data.

NOTE THAT THESE OBJECTS ARE LIKELY TO BE CONVERTED TO S4 OBJECTS, IN THE NEAR FUTURE.\cr
ALSO, NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.

IN PRINCIPLE, YOU SHOULD NOT ACCESS ATTRIBUTES/SLOTS, DIRECTLY.
}
\usage{
pmfuv.dks (x = 1:length (h), h=1, \dots,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default",
    Xlim = c (a, b),
    a = min1 (x), b=Inf)

cdfuv.dks (x = 1:length (h), h=1, \dots,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", tail="lower",
    Xlim = c (a, b),
    a = min1 (x), b=Inf)

qfuv.dks (x = 1:length (h), h=1, \dots,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default",
    Xlim = c (a, b),
    a = min1 (x), b=Inf)
}
\arguments{
\item{x}{
	Integer vector of integer-indexed discrete observations, or bins of such observations.\cr
	(If duplicates, they, along with their frequencies, are aggregated).\cr 
	Also, can be a single-column integer matrix, preferably with a column (variable) name, and optionally with row (bin) names.\cr\cr
	Defaults to a sequence from one to the length of h.}
\item{h}{
	Positive numeric vector of frequencies (or weights), which can be fractional.\cr
	(If scalar, it's recycled to match the length of pre-aggregate x).\cr\cr
	Defaults to one, such that each x value represents a single discrete value.}
\item{bw}{
	Odd positive integer value, giving the bandwidth parameter.\cr
	If missing, an initial bandwidth is computed via the the bandwidth method (see bw.method below), which is subject to the smoothness parameter.}
\item{smoothness}{
	Positive numeric value, giving the relative bandwidth.\cr
	Ignored, if bw is provided.}
\item{kernel}{A (continuous) \link{kernel} object.}
\item{bw.method}{String, the bandwidth selection method.\cr Refer to \link{Bandwidth Selection}.}
\item{tail}{
	String, either "lower" or "upper".\cr
	If lower (the default), lower tail probabilities, P (X <= x), are used.\cr
	If upper, upper tail probabilities, P (X >= x), are used.}
\item{Xlim}{
	In principle, a length-two integer vector, giving the limits of X.\cr
	But a numeric vector is allowed, to support -Inf/Inf.\cr\cr
	The corresponding random variable is regarded as bounded, if either limit is finite.\cr
	In which case, a truncated smoothing algorithm is applied.}
\item{a, b}{
	In principle, integer values.\cr
	This is an alternative way of specifying Xlim, above.\cr\cr
	The min1 function will return one, if one is the minimum x value, otherwise, it will return zero.}
\item{\dots}{Additional arguments not allowed.}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

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
\link{Kernels}

\link{Succinct Constructors}\cr
\link{Continuous Kernel Smoothing}, \link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{is.dks}, \link{ph.printf.phmodel}, \link{ph.plotf.dksuv}

\link{Bandwidth Selection}
}
\examples{
prep.ph.data ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq)
dFh <- cdfuv.dks (traffic.bins, traffic.freq)
dFht <- qfuv.dks (traffic.bins, traffic.freq)

plot (dfh)
plot (dfh, TRUE)
plot (dfh, freq=TRUE)
plot (dFh, freq=TRUE)

dFht (0.5)
}
