\name{53_pairwise_kernel_arrays}
\alias{list.ckernels}
\alias{plot_kernel_array}
\title{Pairwise Kernel Arrays}
\description{Plots of pairwise kernel arrays.}
\usage{
list.ckernels ()

plot_kernel_array (ks = list.ckernels (), \dots, ref.line=TRUE, colors)
}
\arguments{
\item{ks}{List of Kernel objects.}
\item{ref.line}{If true, add a reference line.}
\item{colors}{Optional character vector of colors for each plot.}
\item{\dots}{Ignored.}
}
\details{
The function list.ckernels, simply returns a list of CKernel (continuous kernel) objects.

The plot_kernel_array function plot pairs of kernels, for comparison purposes.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Kernels}\cr

\link{ph.plotf.DKernel}, \link{ph.plotf.CKernel}
}
\examples{
plot_kernel_array ()
}
