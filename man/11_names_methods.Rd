\name{11_names_methods}
\alias{ph.namesf.phmodel}
\title{Names Methods}
\description{Print variable names, taken from the input data.}
\usage{
\method{ph.namesf}{phmodel}(sf, \dots, all=FALSE)
}
\arguments{
\item{sf}{An object extending phmodel, which is most of the objects in this package.}
\item{all}{Logical, if false (the default), only include conditional/random variables.}
\item{\dots}{Ignored.}
}
\value{A character vector of variable names.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
prep.ph.data ()

cfh <- pdfmv.cks (trees2)
names (cfh)
}
