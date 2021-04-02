\name{22_categorical_distributions}
\alias{CAT}
\alias{Categorical Distributions}
\alias{pmfuv.cat}
\alias{cdfuv.cat}
\alias{qfuv.cat}
\alias{pmfc.cat}
\alias{cdfc.cat}
\alias{qfc.cat}
\title{Categorical Models}
\description{
Fit categorical distributions, from data.

NOTE THAT THESE OBJECTS ARE LIKELY TO BE CONVERTED TO S4 OBJECTS, IN THE NEAR FUTURE.\cr
ALSO, NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.

IN PRINCIPLE, YOU SHOULD NOT ACCESS ATTRIBUTES/SLOTS, DIRECTLY.
}
\usage{
#univariate
pmfuv.cat (g, h=1)
cdfuv.cat (g, h=1)
qfuv.cat (g, h=1)

#conditional
pmfc.cat (g, h=1, \dots, conditions, warning=TRUE)
cdfc.cat (g, h=1, \dots, conditions, warning=TRUE)
qfc.cat (g, h=1, \dots, conditions, warning=TRUE)
}
\arguments{
\item{g}{Integer/factor/character vector of groups.\cr Also, can be a named list of such vectors.\cr\cr
	For univariate distributions, the list should only have one vector.\cr
	For conditional distributions, the list needs two or more equal-length vectors.}
\item{h}{
	Optional numeric vector of frequencies (or weights).\cr
	It's length should be one or n, equal the length of the g vectors.}
\item{conditions}{
	An integer vector of category indices, a character vector of category names, or a list which can contain either integers (indices) or strings (names).\cr
	The vector or list can be named (which is preferable) or unnamed.\cr\cr
	If named, then the names are matched against the variable names.\cr
	If unnamed, then the first condition applies to the first variable, and the second condition applies to the second variable, and so on.\cr\cr
	Note that the number of conditions needs to equal the number of variables minus one.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{Additional arguments not allowed.}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

Refer to the vignette for more information.

Default variable names are generated, if the g list is unnamed.\cr
(In conditional models, a warning is generated).
}
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
e.g. If unique g values were 10, 11 and 12, in contrast to 1, 2 and 3.

}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}\cr
These can be used for statistical classification purposes.

\link{is.cat}, \link{ph.printf.phmodel}, \link{ph.plotf.catuv}
}
\examples{
prep.ph.data ()

gfh <- pmfuv.cat (crime.type, n.arrests)
gFht <- qfuv.cat (crime.type, n.arrests)

plot (gfh, freq=TRUE)
ph.mode (gfh)
ph.mode (gfh, level.names=TRUE)

gFht (0.5)
gFht (0.5, category=TRUE)
}
