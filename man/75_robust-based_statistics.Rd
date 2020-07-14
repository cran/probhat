\name{75_robust-based_statistics}
\alias{ph.median}
\alias{ph.quantile}
\alias{iqr}
\title{Robust Statistics}
\description{Compute robust statistics from probability distributions.}
\usage{
ph.median (xf, \dots)
ph.quantile (xf, p, \dots)

iqr (xf, P=0.5, \dots)
}
\arguments{
\item{xf}{A numeric vector, suitable function object, or an object that can be coerced to a numeric vector.\cr Here, a suitable function object is a quantile function.\cr\cr Refer to the references and see also sections.}
\item{P, p}{Numeric vectors, the probabilities.\cr P is the area (probability) between the lower and upper limits.}
\item{\dots}{Other arguments.\cr Refer to the details section.}
}
\details{
If xf is a numeric vector, a \link{qfuv.el} object is created using xf as the main argument.\cr
Any arguments contained within \dots, are passed to the qfuv.el constructor.

If xf is not a quantile function, these functions try to coerce it to a numeric vector, and apply the above.
}
\value{
ph.median returns a single numeric value.

The other functions return a numeric vector.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Empirical-Like Distributions}

\link{probmv}, \link{rng}

\link{ph.mean}, \link{moment}\cr
\link{quartiles}, \link{ntiles}\cr
\link{ph.mode}, \link{ph.modes}
}
\examples{
ph.data.prep ()

cFht <- qfuv.cks (height)

cFht (0.5)
ph.median (cFht)

iqr (cFht)
}
