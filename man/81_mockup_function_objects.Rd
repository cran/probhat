\name{81_mockup_function_objects}
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
\description{Hard-coded functions, representing (runtime) function objects.

DO NO CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.}
\usage{
#discrete kernel smoothing models
#(DKS)
mf.dfh (x, \dots, freq)         #PMF
mf.dFh (q)                    #CDF
mf.dFht (p)                   #QF

#categorical models
#(CAT/gMIXp)
mf.gfh (g, \dots, freq)         #PMF
mf.gFh (q)                    #CDF
mf.gFht (p, \dots, name=FALSE)  #QF

#univariate continuous models
#(CKS/EL/xMIXp, UV/C)
mf.cfh (x)                    #PDF
mf.cFh (q)                    #CDF
mf.cFht (p)                   #QF

#multivariate continuous models
#(CKS, MV/MVC)
mf.cfh.mv (x)                 #PDF
mf.cFh.mv (q)                 #CDF

#chained quantile functions
mf.chFht (p)
}
\arguments{
\item{g, x, q}{
	IN (DKS) MODELS:\cr
	An integer vector, of quantiles.\cr\cr
	IN (CAT/gMIXp) MODELS:\cr
	An integer/factor/character vector, of quantiles.\cr Integers represent category indices.\cr Characters and formatted factors represent category names.\cr\cr
	IN (CKS/EL/xMIXp, UV/C) MODELS:\cr
	A numeric vector, of quantiles.\cr\cr
	IN (CKS, MV/MVC) MODELS:\cr
	A numeric vector or matrix, of quantiles.\cr
	Standard numeric vectors are rbind-ed into single-row matrices.\cr
	Each row represents one evaluation point, and each column represents one variable.}
\item{p}{A numeric vector of probabilities, between zero and one.\cr\cr Except in chained quantile functions, where p should be a numeric vector or matrix.\cr Standard numeric vectors are rbind-ed into single-row matrices.}
\item{freq}{Logical, if true, return frequencies rather than probabilities.\cr Refer to the value section.\cr\cr Default depends on the object.\cr Refer to the corresponding constructors.}
\item{name}{Logical, if true, return category names rather than category indices.\cr Refer to the value section.}
\item{\dots}{.}
}
\details{DO NO CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.

If x, q or p are matrices, then the order of the columns should be the same as the order of the random variables in the model.}
\value{
PMFs return a numeric vector, giving mass (if freq=FALSE) or frequencies (if freq=TRUE).

PDFs return a numeric vector, giving density.

CDFs return a numeric vector, giving cumulative probability (from zero to one).

In (DKS) models, QFs return an integer vector of quantiles.

In (CAT/gMIXp) models, QFs return an integer vector (if name=FALSE) of category indices or a character vector (if name=TRUE) of category names.

In (CKS/EL/xMIXp, UV/C) models, QFs return a numeric vector of quantiles.

Chained quantiles, return a numeric matrix of multivariate quantiles.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}
}

