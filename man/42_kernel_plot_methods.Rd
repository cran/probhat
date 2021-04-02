\name{42_kernel_plot_methods}
\alias{ph.plotf.DKernel}
\alias{ph.plotf.CKernel}
\title{Kernel Plot Methods}
\description{Plot methods for kernel objects.}
\usage{
####################################
#discrete kernels
#(call plot_dpd)
####################################
\method{ph.plotf}{DKernel}(dk, \dots, cdf=FALSE)

####################################
#continuous kernels
#(call plot_cpd)
####################################
\method{ph.plotf}{CKernel}(ck, \dots, cdf=FALSE)
}
\arguments{
\item{dk, ck}{A kernel object.\cr Refer to the references and see also sections.}
\item{cdf}{Logical, if true, plot the CDF.}
\item{\dots}{Other arguments, for plot_dpd or plot_cpd.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Kernels}

\link{plot_kernel_array}

\link{plot_dpd}, \link{plot_cpd}
}
\examples{
dk <- discretized.kernel (7)
plot (dk)

plot (BIWEIGHT.CKERNEL)
}
