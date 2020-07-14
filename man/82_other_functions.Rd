\name{82_other_functions}
\alias{ph.data.prep}
\title{Other Functions}
\description{Emulate an R script, and prepare data for examples.}
\usage{
ph.data.prep (\dots, eval=TRUE, echo=FALSE)
}
\arguments{
\item{eval}{Logical, if true, run the script.}
\item{echo}{Logical, if true, echo the script.}
\item{\dots}{Ignored.}
}
\details{
The ph.data.prep function, which is only designed for tests and examples, runs an R script stored inside the function, itself.
}
\references{Refer to the vignette for an overview, references and better examples.}
\examples{
ph.data.prep (eval=FALSE, echo=TRUE)
}
