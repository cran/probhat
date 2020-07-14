\name{72_random_number_generation}
\alias{rng}
\title{Random Numbers}
\description{Generate random numbers (or synthetic data), univariate or multivariate.}
\usage{
rng (xf, n=1, \dots)
}

\arguments{
\item{xf}{A numeric vector, suitable function object, or an object that can be coerced to a numeric vector.\cr Here, a suitable function object is quantile function, or a chained quantile function.\cr\cr Refer to the references and see also sections.}
\item{n}{Integer, number of random numbers.}
\item{\dots}{Other arguments.\cr\cr Refer to the details section.}
}
\details{
If xf is a numeric vector, a \link{qfuv.el} object is created using xf as the main argument.\cr
Any arguments contained within \dots, are passed to the qfuv.el constructor.

If xf is not a quantile function, these functions try to coerce it to a numeric vector, and apply the above.

Note that the method used for multivariate random number generation is not efficient.
}
\value{
A numeric vector, or numeric matrix.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
ph.data.prep ()

cFht <- qfuv.cks (height)
rng (cFht, 30)

chFht <- chqf.cks (trees)
rng (chFht, 30)

rng (height, 30)
}
