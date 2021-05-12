\name{24_empirical-like_distributions}
\alias{Empirical-Like Distributions}
\alias{cdfuv.el}
\alias{qfuv.el}
\title{Empirical-Like Models}
\description{
Fit empirical-like probability distributions, from data.

NOTE THAT THESE OBJECTS ARE LIKELY TO BE CONVERTED TO S4 OBJECTS, IN THE NEAR FUTURE.\cr
ALSO, NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.

IN PRINCIPLE, YOU SHOULD NOT ACCESS ATTRIBUTES/SLOTS, DIRECTLY.
}
\usage{
cdfuv.el (x, \dots, w)
qfuv.el (x, \dots, w)
}
\arguments{
\item{x}{Numeric vector of data.\cr Also can be a single-column numeric matrix, preferably with a column name.}
\item{w}{Optional numeric vector of weights.}
\item{\dots}{Additional arguments not allowed.}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

Refer to the vignette for more information.}
\value{
Self-referencing function objects.

Refer to \link{Runtime Function Objects}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Categorical Distributions}

\link{is.el}, \link{ph.printf.phmodel}, \link{ph.plotf.eluv}
}
\examples{
prep.ph.data ()

eFht <- qfuv.el (height)
eFht (c (0.25, 0.5, 0.75) )
}
