\name{02_standardized_continuous_kernels}
\alias{uniform.kernel}
\alias{triangular.kernel}
\alias{epanechnikov.kernel}
\alias{truncnorm.kernel}
\alias{biweight.kernel}
\alias{triweight.kernel}
\alias{tricube.kernel}
\alias{bell.spline}
\title{Standardized Continuous Kernels}
\description{Standardized continuous kernels for continuous kernel smoothing.}
\usage{
uniform.kernel ()
triangular.kernel ()
epanechnikov.kernel ()
truncnorm.kernel ()
biweight.kernel ()
triweight.kernel ()
tricube.kernel ()
bell.spline ()
}
\details{
These are constructors for standardized continuous kernel objects.

Currently, these objects are named lists, however, it's possible that they may be changed to S4 objects, in the future.

Each object has $pdf and $cdf components, each of which are functions, giving its PDF and CDF, respectively.

Here, PDFs are symmetric about zero, and have zero-density outside the interval [-1, 1].

The truncnorm constructor creates a truncated normal (or truncated Gaussian) kernel, which is symmetrically truncated, such that the truncated area from the untruncated distribuion, is 0.995 of the untruncated distribution. It's possible that I may allow relatively arbitrary truncation points, in the future.

The bell spline is a novel kernel, constructed from a three-piece quadratic spline, with one control point at (-0.5, 0.5) and another at (0.5, 0.5).

Note that the kernel.array function can be used to plot and compare multiple kernels.
}
\value{
A named list with two component functions.

Each component function has the form:

function (x) = ...

Where x is a numeric vector.
}
\examples{
ckernel = epanechnikov.kernel ()
ckernel$pdf (0)
ckernel$cdf (0)
}
