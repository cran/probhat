\name{33_print_methods}
\alias{print.phmodel}
\title{Print Methods}
\description{Print methods for objects in this package.}
\usage{
\method{print}{phmodel}(x, \dots)
}
\arguments{
\item{x}{An object extending phmodel, which is most of the objects in this package.}
\item{\dots}{Other arguments for the object.summary function.}
}
\details{The print method calls intoo::object.summary.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
ph.data.prep ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq, lower=0)
dfh
}
