\name{32_range_and_sequence_methods}
\alias{min.dpduv}
\alias{max.dpduv}
\alias{min.cpduv}
\alias{max.cpduv}
\alias{range.dpduv}
\alias{range.cpduv}
\alias{seq.dpduv}
\alias{seq.cpduv}
\title{Range and Sequence Functions}
\description{Range and sequence methods for probability distributions.}
\usage{
\method{min}{dpduv}(sf, infv=FALSE, \dots)
\method{max}{dpduv}(sf, infv=FALSE, \dots)
\method{min}{cpduv}(sf, infv=FALSE, \dots)
\method{max}{cpduv}(sf, infv=FALSE, \dots)

\method{range}{dpduv}(sf, infv=FALSE, \dots, freq)
\method{range}{cpduv}(sf, infv=FALSE, \dots)

\method{seq}{dpduv}(sf, infv=FALSE, \dots, midpoints=TRUE, freq)
\method{seq}{cpduv}(sf, infv=FALSE, \dots, n=200)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this refers to a univariate probability distribution.\cr Refer to the references and see also sections.}
\item{infv}{Logical, in function value.\cr Except for quantile functions, where this refers to the probabilities.\cr Refer to the details section.}
\item{midpoints}{Logical, if true, return midpoints.\cr Ignored, except for discrete quantile functions with infv=TRUE.\cr Refer to the details section.}
\item{freq}{Logical, if true, return frequencies.\cr Ignored, except for PMFs with infv=TRUE.}
\item{n}{Integer, the number of points in the resulting sequence.}
\item{\dots}{Ignored.}
}
\details{
By default, the min/max, range and sequence methods apply to range of the random variable.\cr
Often this the range of the observations plus/minus half the bandwidth at each end.

Calling the sequence method on a discrete quantile function, with infv=TRUE:

If midpoints is true, then midpoints of the intervals are returned.\cr
If midpoints is false, then breakpoints, including the outermost values, are returned.

Each interval is defined by one consecutive pair of breakpoints.\cr
Where the breakpoints are (unique) values from the CDF, including zero (at the start) and one (at the end).

In general, these sequences are not equally-spaced.

Calling the sequence method on a continuous quantile function, with infv=TRUE:

Simply returns a returns an equally-spaced sequence between zero and one.
}
\value{
Integer types are returned for discrete probability distributions with infv=FALSE.\cr
Otherwise, numeric types are returned.

The min and max methods return a single integer/numeric value.\cr
The range methods return an length-two integer/numeric vector.

And the seq methods return an integer/numeric vector.\cr
This will be equally-spaced, if infv=FALSE.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
ph.data.prep ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq, lower=0)

seq (dfh)
seq (dfh, TRUE)
}
