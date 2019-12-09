\name{06_empirical_like_distributions}
\alias{cdf.el}
\alias{qf.el}
\title{Empirical-Like Distributions}
\description{Fit (estimated) empirical-like distributions, from data.}
\usage{
cdf.el (x, w=NA)
qf.el (x, w=NA)
}
\details{Refer to the vignette for more information.}
\value{
Self-referencing functions of the form:

function (x) = ...\cr
function (y) = ...

Where x and y are numeric vectors, and y has values in the interval [0, 1].
}
\arguments{
\item{x}{Numeric vector or a single-column numeric matrix, preferably with a column name.}
\item{w}{Optional numeric vector of weights.}
}
\examples{
data.prep ()

Fh.inv = qf.el (Height)
Fh.inv (0.5)
}
