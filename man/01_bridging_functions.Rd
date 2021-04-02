\name{01_bridging_functions}
\alias{ph.namesf}
\alias{ph.printf}
\alias{ph.plotf}
\alias{ph.linesf}
\alias{names.phob}
\alias{print.phob}
\alias{plot.phob}
\alias{lines.phob}
\title{Bridging Functions}
\description{
Generic functions for calling S3 generics with function objects, or objects that contain function objects.\cr
In this package, such objects are often labelled, sf, for suitable function.\cr
(They're also often used to represent probability distributions).

In general, it's easier to call the standard generic.\cr
i.e. It's easier to call print() rather than ph.printf().

However, if using the standard generic, it's best not to name the first argument.
}
\usage{
ph.namesf (\dots)
ph.printf (\dots)
ph.plotf (\dots)
ph.linesf (\dots)

\method{names}{phob}(x, \dots)
\method{print}{phob}(x, \dots)
\method{plot}{phob}(x, \dots)
\method{lines}{phob}(x, \dots)
}
\arguments{
\item{x}{The first argument.}
\item{\dots}{Further arguments, for the S3 method.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\examples{
print (BIWEIGHT.CKERNEL)
ph.printf (BIWEIGHT.CKERNEL)
}
