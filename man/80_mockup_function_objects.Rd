\name{80_mockup_function_objects}
\alias{Mockup Function Objects}
\alias{mf.dfh}
\alias{mf.dFh}
\alias{mf.dFht}
\alias{mf.gfh}
\alias{mf.gFh}
\alias{mf.gFht}
\alias{mf.cfh}
\alias{mf.cFh}
\alias{mf.cFht}
\alias{mf.cfh.mv}
\alias{mf.cFh.mv}
\alias{mf.chFht}
\title{Mockup Function Objects}
\description{
Hard-coded functions, representing (runtime) function objects.

DO NOT CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.\cr
THE RESULTING FUNCTION OBJECTS SHOULD HAVE THE SAME ARGUMENTS AS BELOW.}
\usage{
#categorical models
#(CAT/gMIX)
mf.gfh (g, \dots, freq=FALSE, n)      #PMF
mf.gFh (g, \dots, freq=FALSE, n)      #CDF
mf.gFht (p, \dots, level.names=FALSE) #QF

#other discrete models
#(DKS)
mf.dfh (x, \dots, freq=FALSE, n)      #PMF
mf.dFh (x, \dots, freq=FALSE, n)      #CDF
mf.dFht (p)                         #QF

#univariate continuous models
#(CKS/EL/xMIX, UV/C)
mf.cfh (x)                          #PDF
mf.cFh (x)                          #CDF
mf.cFht (p)                         #QF

#multivariate continuous models
#(CKS, MV/MVC)
mf.cfh.mv (x)                       #PDF
mf.cFh.mv (x)                       #CDF

#chained quantile functions
mf.chFht (p)
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

