\name{rd_7_combinatoric_probability}
\alias{comb.prob}
\title{Combinatoric Probability}
\description{Compute probabilities (including multivariate probabilities) by evaluating CDFs, combinatorically.}
\usage{
comb.prob (cdf.f, a, b)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{cdf.f}{A CDF.}
\item{a}{A vector (or matrix) of lower limits. If a matrix, then each row represents one region and each column represents one variable.}
\item{b}{Same as a, except upper limits.}
}
\examples{
data (trees)
trees = as.matrix (trees)
npcdfmv.f = npcdfmv (trees)

summary (trees)

#approximate first and third quartiles
#marginal probabilities are approximately 0.5 each
a = c (11.05, 72, 19.40)
b = c (15.25, 80, 37.30)

#however, multivariate probability is not necessarily close to 0.5 ^ 3
comb.prob (npcdfmv.f, a, b)
}
