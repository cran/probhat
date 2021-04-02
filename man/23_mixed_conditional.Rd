\name{23_mixed_conditional}
\alias{Conditional Distributions with Mixed Input Types}
\alias{pmfc.gmix}
\alias{cdfc.gmix}
\alias{qfc.gmix}
\alias{pdfc.xmix}
\alias{cdfc.xmix}
\alias{qfc.xmix}
\title{Conditional Distributions with Mixed Categorical-Continuous Input}
\description{
	Fit conditional categorical or continuous distributions with mixed categorical-continuous input.\cr
	(These resemble conditional \link{CAT} and \link{CKS} models).\cr
	And the categorical distributions can be used for statistical classification purposes.

	NOTE THAT THESE OBJECTS ARE LIKELY TO BE CONVERTED TO S4 OBJECTS, IN THE NEAR FUTURE.\cr
	ALSO, NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.

	IN PRINCIPLE, YOU SHOULD NOT ACCESS ATTRIBUTES/SLOTS, DIRECTLY.
}
\usage{
#conditional categorical
pmfc.gmix (g, x, \dots, conditions, warning=TRUE, w)
cdfc.gmix (g, x, \dots, conditions, warning=TRUE, w)
qfc.gmix (g, x, \dots, conditions, warning=TRUE, w)

#conditional continuous
pdfc.xmix (g, x, \dots, conditions, warning=TRUE, w)
cdfc.xmix (g, x, \dots, conditions, warning=TRUE, w)
qfc.xmix (g, x, \dots, conditions, warning=TRUE, w)
}
\arguments{
\item{g}{Integer/factor/character vector of groups.\cr Also, can be a named list of one or more such vectors.}
\item{x}{A numeric vector or a numeric matrix, preferably with column names.\cr The length of x (if standard vector) or the number or rows (if a matrix) should equal the length of the g vectors.}
\item{conditions}{Refer to the conditions arg in categorical and continuous conditional models.\cr
	This is the same, except that the vector or list, needs to be named (unnamed conditions are not allowed), and can include both categorical and continuous variables.\cr\cr
	Note that the number of conditions needs to equal the number of variables minus one.\cr
	(For categorical distributions, there should be one categorical variable left out, and for continuous distributions there should be one continuous variable left out).\cr\cr
	The resulting probability distribution should be the conditional distribution of the variable not included in the conditions.}
\item{w}{Optional numeric vector of weights.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{
	In categorical distributions, further arguments for pdfmv.cks, which is called on the continuous conditioning variables.\cr
	In continuous distributions, further arguments for pdfuv or pdfc, which is called on the continuous conditional variable.}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

Refer to the vignette for more information.

Default variable names are generated, if the g/x list/matrix are unnamed.\cr
(And a warning is generated).

Note that categorical and continuous variables need different names.
}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}

Except:\cr
The constructors for conditional distributions, return NULL, if there's no observations within the conditional window.\cr
(And by default, generate a warning).

Note you can check names (and their order), with names method.\cr
This may be useful for conditional distributions.\cr
(By default, only the conditional variables are returned).
}
\note{
WARNING:

If a conditional categorical distribution is constructed with integer g values, the category indices won't necessarily equal the category names.
e.g. If unique g values were 10, 11 and 12, in contrast to 1, 2 and 3.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{is.cat}, \link{ph.printf.phmodel}, \link{ph.plotf.catuv}
}
\examples{
prep.ph.data ()

fh.gmix <- pmfc.gmix (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )
Fht.gmix <- qfc.gmix (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )

plot (fh.gmix)
ph.mode (fh.gmix)
ph.mode (fh.gmix, level.names=TRUE)

Fht.gmix (0.5)
Fht.gmix (0.5, category=TRUE)
}
