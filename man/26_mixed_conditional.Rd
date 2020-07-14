\name{26_mixed_conditional}
\alias{Conditional Distributions with Mixed Input Types}
\alias{pmfc.gmixp}
\alias{cdfc.gmixp}
\alias{qfc.gmixp}
\alias{pdfc.xmixp}
\alias{cdfc.xmixp}
\alias{qfc.xmixp}
\title{Conditional Distributions with Mixed Input Types}
\description{Fit conditional categorical or continuous distributions with mixed input data.\cr Such categorical distributions can be used for statistical classification purposes.}
\usage{
#conditional categorical
pmfc.gmixp (g, x, \dots, conditions, warning=TRUE, w=NA, freq=FALSE)
cdfc.gmixp (g, x, \dots, conditions, warning=TRUE, w=NA)
qfc.gmixp (g, x, \dots, conditions, warning=TRUE, w=NA)

#conditional continuous
pdfc.xmixp (g, x, \dots, conditions, warning=TRUE, w=NA)
cdfc.xmixp (g, x, \dots, conditions, warning=TRUE, w=NA)
qfc.xmixp (g, x, \dots, conditions, warning=TRUE, w=NA)
}
\arguments{
\item{g}{Integer/factor/character vector of groups.\cr Also, can be a named list of one or more such vectors.}
\item{x}{A numeric vector or a numeric matrix, preferably with column names.\cr The length of x (if standard vector) or the number or rows (if a matrix) should equal the length of the g vectors.}
\item{conditions}{Refer to the conditions arg in categorical and continuous conditional models.\cr
	This is the same, except that the vector or list, needs to be named (unnamed conditions are not allowed), and can include both categorical and continuous variables.\cr\cr
	Note that the number of conditions needs to equal the number of variables minus one.\cr
	(For categorical distributions, there should be one categorical variable left out, and for continuous distributions there should be one continuous variable left out).\cr\cr
	The resulting probability distribution should be the conditional distribution of the variable not included in the conditions.}
\item{freq}{Logical, if true, the resulting function object returns frequencies, by default.}
\item{w}{Optional numeric vector of weights.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{
	In categorical distributions, further arguments for pdfmv.cks, which is called on the continuous conditioning variables.\cr
	In continuous distributions, further arguments for pdfuv or pdfc, which is called on the continuous conditional variable.}
}
\details{Refer to the vignette for more information.

Note that categorical and continuous variables need different names.}
\value{
Self-referencing function objects.

Refer to \link{Mockup Function Objects}

Except:\cr
The constructors for conditional distributions, return NULL, if there's no observations within the conditional window.\cr
(And by default, generate a warning).
}
\note{
WARNING:

If a conditional categorical distribution is constructed with integer g values, the category indices won't necessarily equal the category names.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{is.cat}, \link{print.phmodel}, \link{plot.catuv}
}
\examples{
ph.data.prep ()

fh.gmix <- pmfc.gmixp (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )
Fht.gmix <- qfc.gmixp (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )

ph.mode (fh.gmix)
gmode (fh.gmix)

Fht.gmix (0.5)
Fht.gmix (0.5, category=TRUE)
}
