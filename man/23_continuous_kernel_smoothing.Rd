\name{23_continuous_kernel_smoothing}
\alias{Continuous Kernel Smoothing}
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
\title{Continuous Kernel Smoothing Models}
\description{Fit probability distributions, via continuous kernel smoothing, from data.}
\usage{
#univariate
pdfuv.cks (x, \dots, spline=TRUE,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA)
cdfuv.cks (x, \dots, spline=TRUE,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA)
qfuv.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA)

#multivariate
pdfmv.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    bw, smoothness=1, w=NA)
cdfmv.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    bw, smoothness=1, w=NA)

#conditional
pdfc.cks (x, \dots, spline=TRUE,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA,
    conditions, preserve.range=FALSE,
    warning=TRUE)
cdfc.cks (x, \dots, spline=TRUE,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA,
    conditions, preserve.range=FALSE,
    warning=TRUE)
qfc.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=30, bw, smoothness=1, w=NA,
    conditions, preserve.range=FALSE,
    warning=TRUE)

#multivariate-conditional
pdfmvc.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    bw, smoothness=1, w=NA,
    conditions, preserve.range=FALSE,
    warning=TRUE)
cdfmvc.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    bw, smoothness=1, w=NA,
    conditions, preserve.range=FALSE,
    warning=TRUE)

#other
chqf.cks (x, \dots,
    bw.method="ph.default",
    kernel=biweight.ckernel,
    nc=16, bw, smoothness=1, w=NA)
}
\arguments{
\item{x}{
	IN UNIVARIATE CONSTRUCTORS:\cr
	Numeric vector of data.\cr
	Also, can be a single-column numeric matrix, preferably with a column name.

	IN OTHER CONSTRUCTORS:\cr
	Numeric matrix, preferably, with column names.}
\item{spline}{Logical, if true, use cubic Hermite splines as intermediate models.\cr In general, this should be true.}
\item{bw.method}{String, the bandwidth selection method.\cr Refer to \link{Bandwidth Selection}.}
\item{kernel}{Constructor for a ckernel (continuous kernel) object.}
\item{nc}{Integer, number of control points, in the spline.\cr Ignored, if spline is false.}
\item{bw, smoothness}{Numeric vectors of length one, or length m (the number of variables), the bandwidth and smoothness parameters.\cr
	If bw is missing, the bandwidth is computed for each variable using bw.method (see above) and the smoothness.}
\item{w}{Optional numeric vector of weights.}
\item{conditions}{
	A numeric vector of conditioning values, preferably named.\cr\cr
	If named, then the names are matched against the variable names.\cr
	If unnamed, then the first condition applies to the first variable, and the second condition applies to the second variable, and so on.\cr\cr
	Note that in univariate-conditional distributions, the number of conditions needs to equal the number of variables minus one.}
\item{preserve.range}{Logical.\cr
	If true, the default range used for range/sequence methods and plotting functions, will be the same as the original data.\cr
	If false (the default), the range is based on the data within a conditioning window, which is often smaller.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{Ignored.}
}
\details{Refer to the vignette for more information.}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}

Except:\cr
The constructors for conditional distributions, return NULL, if there's no observations within the conditional window.\cr
(And by default, generate a warning).
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Continuous Kernels}

\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}\cr
These can be used for statistical classification purposes.

\link{is.cks}, \link{print.phmodel}, \link{plot.cksuv}, \link{plot.cksmv}

\link{Bandwidth Selection}
}
\examples{
ph.data.prep ()

cFht <- qfuv.cks (height)
cFht (0.5)
}
