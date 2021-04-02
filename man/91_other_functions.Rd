\name{91_other_functions}
\alias{ph.reflect.demo}
\alias{prep.ph.data}
\alias{min1}
\title{Other Functions}
\description{Other functions, for the package vignette, or to simplify examples.}
\usage{
ph.reflect.demo (both.dirs=FALSE)
prep.ph.data (\dots, eval=TRUE, echo=FALSE)
min1 (x)
}
\arguments{
\item{both.dirs}{Logical, if true, reflect in both directions.}
\item{eval}{Logical, if true, run the script.}
\item{echo}{Logical, if true, echo the script.}
\item{x}{An integer vector.}
\item{\dots}{Ignored.}
}
\details{
The ph.reflect.demo function, is designed specifically for the package vignette.\cr
(And it's output is subject to change).

The prep.ph.data function, which is only designed for tests and examples, runs an R script stored inside the function, itself.

The min1 function will return one, if one is the minimum x value, otherwise, it will return zero.}
\references{Refer to the vignette for an overview, references and better examples.}
\examples{
prep.ph.data (eval=FALSE, echo=TRUE)
}
