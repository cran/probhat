\name{25_empirical-like_distributions}
\alias{Empirical-Like Distributions}
\alias{cdfuv.el}
\alias{qfuv.el}
\title{Empirical-Like Models}
\description{Fit empirical-like probability distributions, from data.}
\usage{
cdfuv.el (x, \dots, w=NA)
qfuv.el (x, \dots, w=NA)
}
\arguments{
\item{x}{Numeric vector of data.\cr Also can be a single-column numeric matrix, preferably with a column name.}
\item{w}{Optional numeric vector of weights.}
\item{\dots}{Ignored.}
}
\details{Refer to the vignette for more information.}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Categorical Distributions}

\link{is.el}, \link{print.phmodel}, \link{plot.eluv}
}
\examples{
ph.data.prep ()

eFht <- qfuv.el (height)
eFht (0.5)
}
