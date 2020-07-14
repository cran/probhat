\name{24_categorical_distributions}
\alias{Categorical Distributions}
\alias{pmfuv.cat}
\alias{cdfuv.cat}
\alias{qfuv.cat}
\alias{pmfc.cat}
\alias{cdfc.cat}
\alias{qfc.cat}
\title{Categorical Models}
\description{Fit categorical distributions, from data.}
\usage{
#univariate
pmfuv.cat (g, h, \dots, freq=FALSE)
cdfuv.cat (g, h, \dots)
qfuv.cat (g, h, \dots)

#conditional
pmfc.cat (g, h, \dots, conditions, warning=TRUE, freq=FALSE)
cdfc.cat (g, h, \dots, conditions, warning=TRUE)
qfc.cat (g, h, \dots, conditions, warning=TRUE)
}
\arguments{
\item{g}{Integer/factor/character vector of groups.\cr Also, can be a named list of such vectors.\cr\cr
	For univariate distributions, the list should only have one vector.\cr
	For conditional distributions, the list needs two or more equal-length vectors.}
\item{h}{Optional numeric vector of frequencies (or weights).\cr The length of n should equal the length of the g vectors.}
\item{conditions}{
	An integer vector of category indices, a character vector of category names, or a list which can contain either integers (indices) or strings (names).\cr
	The vector or list can be named (which is preferable) or unnamed.\cr\cr
	If named, then the names are matched against the variable names.\cr
	If unnamed, then the first condition applies to the first variable, and the second condition applies to the second variable, and so on.\cr\cr
	Note that the number of conditions needs to equal the number of variables minus one.}
\item{freq}{Logical, if true, the resulting function object returns frequencies, by default.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{Ignored.}
}
\details{Refer to the vignette for more information.}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}

Except:\cr
The constructors for conditional distributions, return NULL, if there's no observations within the conditional window.\cr
(And by default, generate a warning).
}
\note{
WARNING:

If a categorical distribution is constructed from integers, the category indices won't necessarily equal the category names.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}\cr
These can be used for statistical classification purposes.

\link{is.cat}, \link{print.phmodel}, \link{plot.catuv}
}
\examples{
ph.data.prep ()

gfh <- pmfuv.cat (crime.type, n.arrests)
gFht <- qfuv.cat (crime.type, n.arrests)

ph.mode (gfh)
gmode (gfh)

gFht (0.5)
gFht (0.5, category=TRUE)
}
