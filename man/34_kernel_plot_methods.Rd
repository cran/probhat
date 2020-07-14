\name{34_kernel_plot_methods}
\alias{plot.dkernel}
\alias{plot.ckernel}
\title{Kernel Plot Methods}
\description{Plot methods for kernel objects.}
\usage{
####################################
#discrete kernels
#(call plot_dpd)
####################################
\method{plot}{dkernel}(x, \dots, cdf=FALSE)

####################################
#continuous kernels
#(call plot_cpd)
####################################
\method{plot}{ckernel}(x, \dots, cdf=FALSE)
}
\arguments{
\item{x}{A kernel object.\cr Refer to the references and see also sections.}
\item{cdf}{Logical, if true, plot the CDF.}
\item{\dots}{Other arguments, for plot_dpd or plot_cpd.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernels}, \link{Continuous Kernels}

\link{kernel.array}

\link{plot_dpd}, \link{plot_cpd}
}
\examples{
plot (biweight.ckernel () )
}
