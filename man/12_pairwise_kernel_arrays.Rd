\name{12_pairwise_kernel_arrays}
\alias{ckernels}
\alias{kernel.array}
\title{Pairwise Kernel Arrays}
\description{Plots of pairwise kernel arrays.}
\usage{
ckernels ()

kernel.array (which = ckernels (), reference.line=TRUE, colors)
}
\details{
The function ckernels, simply returns a character vector with the names of some (not all) of the continuous kernels in this package.

The kernel.array function plot pairs of kernels, for comparison purposes.
}
\arguments{
\item{which}{Character vector giving the names of the (continuous) kernel constructors.}
\item{reference.line}{If true, add a reference line.}
\item{colors}{Colors for each plot.}
}
\examples{
kernel.array ()
}
