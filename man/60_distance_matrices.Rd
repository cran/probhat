\name{60_distance_matrices}
\alias{psv}
\alias{pdist}
\alias{ph4.rdist}
\alias{ph4.pcomp2}
\title{Distance Matrices}
\description{
Compute the distance (or dissimilarity) between pairs of density functions.

The main function, pdist maps a list (of density functions) to a matrix.\cr
The rdist function, converts the above matrix to a two-column sorted data.frame.

Note that these methods are relatively new.\cr
I can not give any guarantee of optimality.\cr
(This applies to the whole package, but in particular to these functions).

Also, the distance values should be not be interpreted, other than for ranking purposes.
}
\usage{
psv (sf)

pdist (sf, \dots, sqrt.mse=TRUE)

ph4.rdist (d, n)

ph4.pcomp2 (fh, gh, sqrt.mse=TRUE, aggregate=TRUE,
    dfh = psv (fh), dgh = psv (gh) )
}
\arguments{
\item{fh, gh}{
	A density function.\cr
	Refer to \link{pdfuv.cks} and \link{pdfmv.cks}.}
\item{sf}{
	In psv, a density function.\cr
	In pdist, a list of density functions.\cr
	Optionally, they can be produced by \link{ph4.pdfuv.gset.cks} and \link{ph4.pdfmv.gset.cks}.}
\item{aggregate}{If true (the default), return the average of the two one-sided distances.}
\item{sqrt.mse}{If true (the default), the square root of the MSE is used, otherwise, the MSE is used.}
\item{dfh, dgh}{
	Numeric vectors, of self-evaluated densities, see details.\cr
	Note that no validation is done on these arguments.}
\item{d}{A distance matrix, as returned by pdist.}
\item{n}{
	Integer, the closest n pairs.\cr
	If missing, all are returned.}
\item{\dots}{Ignored.}
}
\details{
Here, self-evaluated density values are computed by evaluating a density function at it's own data.\cr
(In contrast to arbitrary evaluation points).

And cross-evaluated densities are values are computed by evaluating a density function at another density function's data.

The psv function computes self-evaluated densities.\cr
The ph.pcomp2 function computes one (or two) distances.\cr
With a single distance between the average of two one-sided distances.

If dfh-xf is the self-evaluated density of fh, and dgh-xf is the density of gh evaluated at fh's data:\cr
(i.e. Two density functions are evaluated at the same points, which are the data points of the first density function).

Then a one-sided distance can be computed as the mean squared error (MSE), or it's square root.\cr
Where the MSE is the sum of the squared differences between dfh-xf and dgh-xf, over n.

The second one-sided distance is the same, except that f and g are reversed.
}
\value{
psv returns a numeric vector.

pdist returns a numeric square matrix.

By default, ph4.pcomp2 returns a single value.\cr
(If aggregate is false, then it returns a pair of values).
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{pwith}, \link{probmv}

\link{pdfuv.cks}, \link{pdfmv.cks}\cr
\link{ph4.pdfuv.gset.cks}, \link{ph4.pdfmv.gset.cks}
}
\examples{
prep.ph.data ()

gs1a <- ph4.pdfuv.gset.cks (species, cbind (sepal.length) )
gs1b <- ph4.pdfuv.gset.cks (species, cbind (sepal.width) )
gs2 <- ph4.pdfmv.gset.cks (species, cbind (sepal.length, sepal.width) )

d1a <- pdist (gs1a)
d1b <- pdist (gs1b)
d2 <- pdist (gs2)

#print out distance matrix
#(for bivariate models)
d2

#print out distances, ranked/sorted
ph4.rdist (d1a)
ph4.rdist (d1b)
ph4.rdist (d2)
}
