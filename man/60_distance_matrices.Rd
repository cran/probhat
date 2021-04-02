\name{60_distance_matrices}
\alias{psv}
\alias{pdist}
\alias{ph4.pcomp2}
\title{Distance Matrices}
\description{
Compute the distance (or dissimilarity) between pairs of density functions.

The main function, pdist maps a list (of density functions) to a matrix.
}
\usage{
psv (fh)

pdist (fhs, \dots, sqrt.mse=TRUE)

ph4.pcomp2 (fh, gh, sqrt.mse=TRUE, aggregate=TRUE,
    dfh = psv (fh), dgh = psv (gh) )
}
\arguments{
\item{fh, gh}{
	A density function.\cr
	Refer to \link{pdfuv.cks} and \link{pdfmv.cks}.}
\item{fhs}{
	A list of density functions.\cr
	Optionally, they can be produced by \link{pdfuv.gset.cks} and \link{pdfmv.gset.cks}.}
\item{aggregate}{If true (the default), return the average of the two one-sided distances.}
\item{sqrt.mse}{If true (the default), the square root of the MSE is used, otherwise, the MSE is used.}
\item{dfh, dgh}{
	Numeric vectors, of self-evaluated densities, see details.\cr
	Note that no validation is done on these arguments.}
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
\link{pdfuv.gset.cks}, \link{pdfmv.gset.cks}
}
\examples{
prep.ph.data ()

gs <- pdfuv.gset.cks (species, sepal.length)
pdist (gs)
}
