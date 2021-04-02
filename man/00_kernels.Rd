\name{00_kernels}
\alias{pdf-class}
\alias{cpd-class}
\alias{phob-class}
\alias{cpduv-class}
\alias{dcdfuv-class}
\alias{pdfuv-class}
\alias{dcdf-class}
\alias{ccdfuv-class}
\alias{ccdf-class}
\alias{dpd-class}
\alias{pmf-class}
\alias{phpd-class}
\alias{dpduv-class}
\alias{pmfuv-class}
\alias{Kernels}
\alias{Kernel-class}
\alias{DKernel-class}
\alias{CKernel-class}
\alias{Uniform.CKernel-class}
\alias{Triangular.CKernel-class}
\alias{Biweight.CKernel-class}
\alias{Triweight.CKernel-class}
\alias{Tricube.CKernel-class}
\alias{TrGaussian.CKernel-class}
\alias{Epanechnikov.CKernel-class}
\alias{Bell.Spline-class}
\alias{discretized.kernel}
\alias{UNIFORM.CKERNEL}
\alias{TRIANGULAR.CKERNEL}
\alias{EPANECHNIKOV.CKERNEL}
\alias{TRGAUSSIAN.CKERNEL}
\alias{BIWEIGHT.CKERNEL}
\alias{TRIWEIGHT.CKERNEL}
\alias{TRICUBE.CKERNEL}
\alias{BELL.SPLINE}
\title{Kernel Objects}
\description{
Discrete and continuous kernel objects.

NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.
}
\usage{
discretized.kernel (n, ck=BIWEIGHT.CKERNEL, \dots, xlim)

UNIFORM.CKERNEL
TRIANGULAR.CKERNEL
EPANECHNIKOV.CKERNEL
TRGAUSSIAN.CKERNEL
BIWEIGHT.CKERNEL
TRIWEIGHT.CKERNEL
TRICUBE.CKERNEL

BELL.SPLINE
}
\arguments{
\item{n}{
	Integer, number of bins.\cr
	Needs to be positive and odd.\cr
	Ignored, if xlim provided.}
\item{ck}{A continuous kernel object.}
\item{xlim}{A length two ascending integer vector.}
\item{\dots}{Ignored.}
}
\details{
Kernel objects are S4 objects with two slots representing the corresponding PMF/PDF and CDF.

Continuous kernels are predefined constants.\cr
Discrete kernels are constructed by using the discretized.kernel() function to discretize a predefined continuous kernel.

Currently, constructors for both \link{DKS} and \link{CKS} objects take a continuous kernel object.\cr
(Where the DKS constructors discretize it, internally). 

Here, PDFs are symmetric about zero, and have positive density over the interval (-1, 1).

Currently, the truncated Gaussian kernel is symmetrically truncated (then transformed), such that the area from the untruncated distribution is 0.995. The bell spline is a novel kernel, constructed from a three-piece quadratic spline, with knots at -0.5 and 0.5.

Note that the \link{plot_kernel_array} function can be used to plot and compare multiple kernels.
}
\value{A Kernel object.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{ph.plotf.DKernel}, \link{ph.plotf.CKernel}\cr
\link{plot_kernel_array}

\link{DKS} and \link{CKS} Models\cr
\link{Conditional Distributions with Mixed Input Types}
}
\examples{
dk <- discretized.kernel (7)
plot (dk)

plot (BIWEIGHT.CKERNEL)
}
