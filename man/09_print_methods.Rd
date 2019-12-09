\name{09_print_methods}
\alias{print.phmodel}
\title{Print Methods}
\description{Print methods for objects in this package.}
\usage{
#univariate
\method{print}{phmodel}(x, \dots)
}
\details{
The print method calls intoo::object.summary.
}
\arguments{
\item{x}{An object extending phmodel, which is most of the objects in this package.}
\item{\dots}{Other arguments for the object.summary function.}
}
\examples{
data.prep ()

fh = pmfuv.dks (traffic.x, traffic.h, bw=25)
fh
}
