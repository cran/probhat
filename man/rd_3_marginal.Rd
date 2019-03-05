\name{rd_3_marginal}
\alias{marginal}
\title{Marginal Distributions}
\description{Constructs a marginal object, which is a list of univariate distributions.}
\usage{
marginal (npmv)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{npmv}{An nppdfmv or npcdfmv object.}
}
\examples{
data (trees)
trees = as.matrix (trees)
nppdfmv.f = nppdfmv (trees)
muv = marginal (nppdfmv.f)
muv [[1]]
}
