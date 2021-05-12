\name{80_runtime_function_objects}
\alias{Runtime Function Objects}
\alias{ph.dfh.rtf}
\alias{ph.dFh.rtf}
\alias{ph.dFht.rtf}
\alias{ph.gfh.rtf}
\alias{ph.gFh.rtf}
\alias{ph.gFht.rtf}
\alias{ph.cfh.rtf}
\alias{ph.cFh.rtf}
\alias{ph.cFht.rtf}
\alias{ph.cfh.rtf.mv}
\alias{ph.cFh.rtf.mv}
\alias{ph.chFht.rtf}
\title{Runtime Function Objects}
\description{
Hard-coded functions, representing (runtime) function objects.

DO NOT CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.\cr
THE RESULTING FUNCTION OBJECTS SHOULD HAVE THE SAME ARGUMENTS AS BELOW.}
\usage{
#categorical models
#(CAT/gMIX)
ph.gfh.rtf (g, \dots, freq=FALSE, n)      #PMF
ph.gFh.rtf (g, \dots, freq=FALSE, n)      #CDF
ph.gFht.rtf (p, \dots, level.names=FALSE) #QF

#other discrete models
#(DKS)
ph.dfh.rtf (x, \dots, freq=FALSE, n)      #PMF
ph.dFh.rtf (x, \dots, freq=FALSE, n)      #CDF
ph.dFht.rtf (p)                         #QF

#univariate continuous models
#(CKS/EL/xMIX, UV/C)
ph.cfh.rtf (x)                          #PDF
ph.cFh.rtf (x)                          #CDF
ph.cFht.rtf (p)                         #QF

#multivariate continuous models
#(CKS, MV/MVC)
ph.cfh.rtf.mv (x)                       #PDF
ph.cFh.rtf.mv (x)                       #CDF

#chained quantile functions
ph.chFht.rtf (p)
}
\arguments{
\item{g, x}{
	IN (DKS) MODELS:\cr
	An integer vector, of quantiles.\cr\cr
	IN (CAT/gMIX) MODELS:\cr
	An integer/factor/character vector, of quantiles.\cr Integers represent category indices.\cr Characters and formatted factors represent category names.\cr\cr
	IN (CKS/EL/xMIX, UV/C) MODELS:\cr
	A numeric vector, of quantiles.\cr\cr
	IN (CKS, MV/MVC) MODELS:\cr
	A numeric vector or matrix, of quantiles.\cr
	Standard numeric vectors are rbind-ed into single-row matrices.\cr
	Each row represents one evaluation point, and each column represents one variable.}
\item{p}{A numeric vector of probabilities, between zero and one.\cr\cr Except in chained quantile functions, where p should be a numeric vector or matrix.\cr Standard numeric vectors are rbind-ed into single-row matrices.}
\item{freq}{Logical, if true, return frequencies rather than probabilities.}
\item{level.names}{Logical, if true, return category names rather than category indices.\cr Refer to the value section.}
\item{n}{
	Sample size.\cr
	Ignored, unless freq is true.}
\item{\dots}{Ignored.}
}
\details{
DO NO CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.

If x or p are matrices, then the order of the columns should be the same as the order of the random variables in the model.

Also, if a bounded interval was used, then the quantiles need to be within the limits.
}
\value{
By default:\cr
PMFs return a numeric vector, giving probability mass.\cr
PDFs return a numeric vector, giving probability density.\cr
CDFs return a numeric vector, giving cumulative probability (from zero to one).\cr
Discrete quantile functions, return an integer vector of quantiles.\cr
Continuous quantile functions, return a numeric vector of quantiles.

In the discrete case, setting freq=TRUE, scales the values to match the n, the sample size.\cr
By default, n is the number of observation used, or the sum of the unscaled weights/frequencies.

In categorical quantile functions, setting level.names=TRUE, returns a character vector of levels names.

Chained quantile functions, return a numeric matrix of multivariate quantiles.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}
}

