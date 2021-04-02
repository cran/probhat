\name{73_mode_estimation}
\alias{ph.mode}
\alias{ph.modes}
\title{Mode Estimation}
\description{Compute the mode or modes, from probability distributions.}
\usage{
ph.mode (sf, infv=FALSE, ..., level.names=FALSE, freq, n)
ph.modes (sf, infv=FALSE)
}
\arguments{
\item{sf}{A suitable function object.\cr For ph.mode, this is a univariate PMF or spline-based PDF.\cr For ph.modes, a spline-based PDF only.\cr\cr Refer to the references and see also sections.}
\item{infv}{Logical, if true, return the value of the PMF/PDF at the mode(s).}
\item{level.names}{Logical, if false, return the category index, if true, return the category name.\cr Ignored, except for categorical PMFs with infv=FALSE.}
\item{freq}{Logical, if true, return frequencies.\cr Ignored, except for PMFs with infv=TRUE.}
\item{n}{
	Integer, the sample size.\cr
	Ignored, unless both infv and freq are true.}
\item{\dots}{Ignored.}
}
\value{
ph.mode returns a single integer or numeric value.\cr
(Except for categorical PMFs with infv=FALSE and level.names=TRUE, which return a string).

ph.modes returns a numeric vector.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing},
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{probmv}, \link{rng}

\link{ph.mean}, \link{moment}\cr
\link{quartiles}, \link{ntiles}, \link{ph.median}, \link{ph.quantile}
}
\examples{
prep.ph.data ()

gfh <- pmfuv.cat (crime.type, n.arrests)
cfh <- pdfuv.cks (height, smoothness=0.5)

ph.mode (gfh)
ph.mode (gfh, TRUE)
ph.mode (gfh, TRUE, freq=TRUE)
ph.mode (gfh, category=TRUE)

ph.mode (cfh)
ph.modes (cfh)
}
