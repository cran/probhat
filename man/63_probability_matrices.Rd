\name{63_probability_matrices}
\alias{BVMatrix-class}
\alias{ph4.BVMatrix}
\title{Bivariate Density/etc Matrices}
\description{
Convenience function to construct an S4 object representing a bivariate density/etc matrix from a bivariate PDF/CDF.

Unlike other objects in the package, you're welcome to access the slots.
}
\usage{
ph4.BVMatrix (sf, xlim, ylim, \dots, n=10)
}
\arguments{
\item{sf}{
	A suitable function object.\cr
	(Here, this refers to a bivariate PDF/CDF).}
\item{xlim, ylim}{Numeric length two vectors, giving the evaluation ranges.}
\item{n}{Integer vector of length one or two, giving the number of points, in each x and y direction.}
\item{\dots}{Ignored.}
}
\value{
An S4 object with three slots, fv (a matrix), x (x evaluation points) and y (y evaluation points).

Each row in fv, corresponds to one x value, and each column to one y value.
}
\references{Refer to the vignette for an overview, references and better examples.}
\examples{
prep.ph.data ()

cfh2 <- pdfmv.cks (trees2 [,-2])
ph4.BVMatrix (cfh2)
}
