\name{71_quantile-based_summaries}
\alias{ntile.names}
\alias{quartiles}
\alias{deciles}
\alias{ntiles}
\title{Order-Based Summary Statistics}
\description{Compute named sequential quantiles, primarily for producing summary-style output.}
\usage{
ntile.names (n, symbol="q", \dots, emph = n / 2)

quartiles (xf, col=FALSE, \dots,
    prob=FALSE, names = ntile.names (4, "Q", emph=emph), emph=2)
deciles (xf, col=FALSE, \dots,
    prob=FALSE, names = ntile.names (10, "D", emph=emph), emph=5)

ntiles (n, xf, col=FALSE, \dots, prob=FALSE, names)
}
\arguments{
\item{n}{Integer, the number of sequential quantiles.}
\item{symbol}{String, letter/symbol for the quantile names.}
\item{xf}{A numeric vector, suitable function object, or an object that can be coerced to a numeric vector.\cr Here, a suitable function object is a quantile function.\cr\cr Refer to the references and see also sections.}
\item{col}{Logical, if true, return a single-column matrix.}
\item{prob}{Logical, if false (the default), name quantiles by index/names, if true, name quantiles by probability.}
\item{names}{Character vector, giving the names.\cr Ignored, if rank is false.}
\item{emph}{In principle, an integer vector in 1:(n-1), which quantiles to emphasize.\cr Can also be a numeric vector, but the floor/ceiling values are used.}
\item{\dots}{Other arguments. \cr Refer to details.}
}
\details{
If xf is a numeric vector, a \link{qfuv.el} object is created using xf as the main argument.\cr
Any arguments contained within \dots, are passed to the qfuv.el constructor.

If xf is not a quantile function, these functions try to coerce it to a numeric vector, and apply the above.
}
\value{
ntiles.names returns a character vector.

The other functions return a named numeric vector (if col=FALSE), or a named single-column matrix (if col=TRUE).
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Empirical-Like Distributions}

\link{probmv}, \link{rng}

\link{ph.mean}, \link{moment}\cr
\link{ph.median}, \link{ph.quantile}\cr
\link{ph.mode}, \link{ph.modes}
}
\examples{
prep.ph.data ()

cFht <- qfuv.cks (height)
quartiles (cFht)
quartiles (cFht, prob=TRUE)

quartiles (height)
}
