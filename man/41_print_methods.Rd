\name{41_print_methods}
\alias{ph.printf.Kernel}
\alias{ph.printf.phmodel}
\alias{ph.printf.dset}
\title{Print Methods}
\description{
Print methods for objects in this package.

Note that this function is likely to change in the next update.
}
\usage{
\method{ph.printf}{Kernel}(k, \dots)
\method{ph.printf}{phmodel}(sf, \dots)
\method{ph.printf}{dset}(vf, \dots)
}
\arguments{
\item{k}{A kernel object.}
\item{sf}{An object extending phmodel, which is most of the objects in this package.}
\item{vf}{A list of probability distributions.}
\item{\dots}{Ignored.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
prep.ph.data ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq)

#autoprinting
dfh
}
