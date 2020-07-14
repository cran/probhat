\name{12_continuous_kernels}
\alias{Continuous Kernels}
\alias{uniform.ckernel}
\alias{triangular.ckernel}
\alias{epanechnikov.ckernel}
\alias{truncnorm.ckernel}
\alias{biweight.ckernel}
\alias{triweight.ckernel}
\alias{tricube.ckernel}
\alias{bell.spline}
\title{Continuous Kernels}
\description{Continuous kernels (over the interval, negative one to positive one), for continuous kernel smoothing.

NOTE THAT THESE MAY BE CHANGED.
}
\usage{
uniform.ckernel ()
triangular.ckernel ()
epanechnikov.ckernel ()
truncnorm.ckernel ()
biweight.ckernel ()
triweight.ckernel ()
tricube.ckernel ()

bell.spline ()
}
\details{
NOTE THAT THESE MAY BE CHANGED.\cr
AND DISCRETE KERNELS MAY BE REMOVED.

These are constructors for continuous kernel objects.

Currently, these objects are named lists (with pdf and cdf components), however, it's possible that they may be changed, in the future.

Here, PDFs are symmetric about zero, and have positive density over the interval (-1, 1).

The truncnorm constructor creates a truncated normal (or truncated Gaussian) kernel, which is symmetrically truncated, such that the area from the untruncated distribution, within the truncation points, is 0.995 of that distribution. It's possible that I may allow relatively arbitrary truncation points, in the future.

The bell spline is a novel kernel, constructed from a three-piece quadratic spline, with knots at -0.5 and 0.5.

These objects are used inside continuous kernel smoothing.\cr
(In which case, you provide the constructor, not the object).

Note that the kernel.array function can be used to plot and compare multiple kernels.
}
\value{A ckernel object.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernels}

\link{plot.ckernel}\cr
\link{kernel.array}

\link{Continuous Kernel Smoothing}
}
\examples{
k <- biweight.ckernel ()
plot (k)

k$pdf (0)
}
