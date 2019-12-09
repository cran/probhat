\name{14_range_and_sequence_functions}
\alias{range.pmfuv}
\alias{range.pdfuv}
\alias{range.cdfuv}
\alias{range.qfuv}
\alias{seq.pmfuv}
\alias{seq.pdfuv}
\alias{seq.cdfuv}
\alias{seq.qfuv}
\title{Range and Sequence Functions}
\description{Range and sequence methods for probability distributions.}
\usage{
\method{range}{pmfuv}(f, \dots)
\method{range}{pdfuv}(f, \dots)
\method{range}{cdfuv}(F, \dots)
\method{range}{qfuv}(F.inv, \dots)

\method{seq}{pmfuv}(f, \dots)
\method{seq}{pdfuv}(f, n, \dots)
\method{seq}{cdfuv}(F, n, \dots)
\method{seq}{qfuv}(F.inv, n, \dots)
}
\details{
These methods apply to both univariate and conditional distributions.\cr

The range relates to the underlying random variable, not the range of the function.

In the context of unbounded kernel smoothing, the range of the random variable is the range of the observed values plus/minus half the bandwidth (or half the bandwidth minus one, in the discrete case).
}
\arguments{
\item{f}{A suitable univariate or conditional distribution.}
\item{F}{A suitable univariate or conditional distribution.}
\item{F.inv}{A suitable univariate or conditional distribution.}
\item{n}{For discrete CDFs, ignored. Otherwise, the number of points in the sequence.}
\item{\dots}{Ignored.}
}
\examples{
data.prep ()

fh = pmfuv.dks (traffic.x, traffic.h, bw=25)
seq (fh)
}
