\name{16_random_numbers}
\alias{ph.rng}
\title{Random Numbers}
\description{Generate random numbers, univariate or multivariate.}
\usage{
ph.rng (F.inv, n)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{F.inv}{A quantile (or chained quantile) function.}
\item{n}{Number of random numbers.}
}
\examples{
data.prep ()

Fh.inv= chqf.cks (trees)
ph.rng (Fh.inv, 10)
}
