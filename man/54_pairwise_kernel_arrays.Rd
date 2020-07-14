\name{54_pairwise_kernel_arrays}
\alias{list.ckernels}
\alias{kernel.array}
\title{Pairwise Kernel Arrays}
\description{Plots of pairwise kernel arrays.}
\usage{
list.ckernels ()

kernel.array (ks = list.ckernels (), \dots, ref.line=TRUE, colors)
}
\arguments{
\item{ks}{List of kernel objects.}
\item{ref.line}{If true, add a reference line.}
\item{colors}{Optional character vector of colors for each plot.}
\item{\dots}{Ignored.}
}
\details{
The function list.ckernels, simply returns a list of ckernel (continuous kernel) objects.

The kernel.array function plot pairs of kernels, for comparison purposes.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernels}, \link{Continuous Kernels}\cr

\link{plot.dkernel}, \link{plot.ckernel}
}
\examples{
kernel.array ()
}
