\name{04_continuous_kernel_smoothing}
\alias{pdfuv.cks}
\alias{cdfuv.cks}
\alias{qfuv.cks}
\alias{pdfmv.cks}
\alias{cdfmv.cks}
\alias{pdfc.cks}
\alias{cdfc.cks}
\alias{qfc.cks}
\alias{pdfmvc.cks}
\alias{cdfmvc.cks}
\alias{chqf.cks}
\title{Models via Continuous Kernel Smoothing}
\description{Fit (estimated) probability distributions, with continuous kernel smoothing, from data.}
\usage{
#univariate
pdfuv.cks (x, spline=TRUE, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
cdfuv.cks (x, spline=TRUE, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
qfuv.cks (x, spline.only, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)

#multivariate
pdfmv.cks (x,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
cdfmv.cks (x,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)

#conditional
pdfc.cks (x, spline=TRUE, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
    conditions, preserve.range=FALSE)
cdfc.cks (x, spline=TRUE, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
    conditions, preserve.range=FALSE)
qfc.cks (x, spline.only, nc=30,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
    conditions, preserve.range=FALSE)

#multivariate-conditional
pdfmvc.cks (x,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
    conditions, preserve.range=FALSE)
cdfmvc.cks (x,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
    conditions, preserve.range=FALSE)

#other
chqf.cks (x, spline.only, nc=16,
    kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
}
\details{
Refer to the vignette for more information.

Expanding on the vignette, conditional models only retain observations within a conditioning window. By default (with preserve.range=FALSE), the range of the random variable is computed from the remaining observations, otherwise, the range is computed from all the observations. The range of the random variable determines the values returned by the range and sequence methods, which in turn, determines the default axes in plots.

If using conditional distributions as an exploratory tool, it may be desirable to set preserve.range to TRUE, because it's likely to be easier to compare different probability distributions.
}
\value{
Self-referencing functions of the form:

function (x) = ...\cr
function (y) = ...

Where, in the univariate case, x and y are numeric vectors, and y has values in the interval [0, 1].

And where, in the multivariate and chained cases, x and y can be either a numeric vector or a numeric matrix. Numeric vectors are turned into numeric matrices by rbind-ing them. And the number of columns must equal the number of random variables.
}
\arguments{
\item{x}{For univariate distributions, x should be a numeric vector or a single-column numeric matrix, preferably with a column name. For multivariate and conditional distributions, it should be a numeric matrix, preferably, with column names.}
\item{spline}{If true, use cubic Hermite splines as intermediate models.}
\item{spline.only}{Ignored.}
\item{nc}{Number of control points, in the spline.}
\item{kernel}{Constructor for a continuous kernel object, or the kernel object itself.}
\item{bw}{The bandwidth parameter of length one or m (the number of variables). If supplied, the smoothness parameter is recomputed, otherwise, the bandwidth is computed from the smoothness.}
\item{smoothness}{A smoothness parameter (equal to the ratio between the bandwidth and the range of the observations) of length one or m (the number of variables).}
\item{w}{Optional numeric vector of weights.}
\item{conditions}{A numeric vector of conditioning values. If it's a named vector, then the names are matched against the names of the variables, otherwise, the first conditioning value applies to the first variable, the second conditioning value applies to the second variable, and so on. For univariate-conditional distributions, the number of conditions should equal m - 1, and for multivariate-conditional distributions, the number of conditions should be in the interval [0, m - 1].}
\item{preserve.range}{Refer to details.}
}
\examples{
data.prep ()

fh = pdfuv.cks (Height)
fh (22)
}
