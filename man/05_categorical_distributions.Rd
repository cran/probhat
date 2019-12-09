\name{05_categorical_distributions}
\alias{pmfuv.cat}
\alias{cdfuv.cat}
\alias{qfuv.cat}
\alias{pmfc.cat.cks}
\title{Categorical distributions}
\description{Fit (estimated) categorical distributions, from data.}
\usage{
#unconditional
pmfuv.cat (x, sequence, w=NA, levels)
cdfuv.cat (x, sequence, w=NA, levels)
qfuv.cat (x, sequence, w=NA, levels)

#conditional
#(on a continuous variable)
pmfc.cat.cks (x, y, sequence,
    kernel=biweight.kernel, bw, smoothness=0.65, w=NA,
    levels, at)
}
\details{Refer to the vignette for more information.
}
\value{
Self-referencing functions of the form:

function (x) = ...\cr
function (y) = ...

Where x is an integer-valued vector, and y is a numeric vector, with values in the interval [0, 1].
}
\arguments{
\item{x}{Integer-valued or character vector, or a single-column integer-valued or character matrix, preferably with a column name.}
\item{y}{Same as x in pdfuv.cks.}
\item{sequence}{Ignored.}
\item{w}{Optional numeric vector of weights.}
\item{levels}{If x is an integer, optional character vector for the levels.}
\item{kernel}{Same as pdfuv.cks.}
\item{bw}{Same as pdfuv.cks.}
\item{smoothness}{Same as pdfuv.cks.}
\item{at}{Conditioning value for the continuous variable.}
}
\examples{
data.prep ()

fh = pmfuv.cat (region)
fh (2)
fh ("Northeast")
}
