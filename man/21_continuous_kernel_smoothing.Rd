\name{21_continuous_kernel_smoothing}
\alias{CKS}
\alias{Continuous Kernel Smoothing}
\alias{pdfuv.cks}
\alias{cdfuv.cks}
\alias{qfuv.cks}
\alias{pdfmv.cks}
\alias{cdfmv.cks}
\alias{pdfc.cks}
\alias{cdfc.cks}
\alias{qfc.cks}
\alias{pdfmvc.cks}
\alias{cdfmvc.cks}
\alias{chqf.cks}
\title{Continuous Kernel Smoothing Models}
\description{
Fit probability distributions, via continuous kernel smoothing, from data.

NOTE THAT THESE OBJECTS ARE LIKELY TO BE CONVERTED TO S4 OBJECTS, IN THE NEAR FUTURE.\cr
ALSO, NOTE THAT THEIR INTERNAL STRUCTURE (THAT IS, THEIR ATTRIBUTES/SLOTS), IS SUBJECT TO CHANGE.

IN PRINCIPLE, YOU SHOULD NOT ACCESS ATTRIBUTES/SLOTS, DIRECTLY.
}
\usage{
#univariate
pdfuv.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    spline=TRUE, bw.method="ph.default", nc=30,
    trtype="local",
    Xlim = cbind (a, b), a=-Inf, b=Inf)

cdfuv.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    spline=TRUE, bw.method="ph.default", nc=30, tail="lower",
    trtype="local",
    Xlim = cbind (a, b), a=-Inf, b=Inf)

qfuv.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", nc=30,
    trtype="local",
    Xlim = cbind (a, b), a=-Inf, b=Inf)

#multivariate
pdfmv.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default",
    Xlim = cbind (a, b), a=-Inf, b=Inf)

cdfmv.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", tail="lower",
    Xlim = cbind (a, b), a=-Inf, b=Inf)

#conditional
pdfc.cks (x, \dots, conditions, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    spline=TRUE, bw.method="ph.default", nc=30,
    Xlim = cbind (a, b), a=-Inf, b=Inf,
    preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE,
    warning=TRUE)

cdfc.cks (x, \dots, conditions, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    spline=TRUE, bw.method="ph.default", nc=30, tail="lower",
    Xlim = cbind (a, b), a=-Inf, b=Inf,
    preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE,
    warning=TRUE)

qfc.cks (x, \dots, conditions, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", nc=30,
    Xlim = cbind (a, b), a=-Inf, b=Inf,
    preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE,
    warning=TRUE)

pdfmvc.cks (x, \dots, conditions, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default",
    Xlim = cbind (a, b), a=-Inf, b=Inf,
    preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE,
    warning=TRUE)

#multivariate-conditional
cdfmvc.cks (x, \dots, conditions, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", tail="lower",
    Xlim = cbind (a, b), a=-Inf, b=Inf,
    preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE,
    warning=TRUE)

#other
chqf.cks (x, \dots, w,
    bw, smoothness=1,
    kernel=BIWEIGHT.CKERNEL,
    bw.method="ph.default", nc=16)
}
\arguments{
\item{x}{
	IN UNIVARIATE CONSTRUCTORS:\cr
	Numeric vector of data.\cr
	Also, can be a single-column numeric matrix, preferably with a column name.

	IN OTHER CONSTRUCTORS:\cr
	Numeric matrix, preferably, with column names.}
\item{conditions}{
	SIMPLE USE:\cr
	(SINGLE SET OF CONDITIONS)\cr
	A numeric vector of conditioning values.\cr\cr
	If named, then the names are matched against the variable names.\cr
	If unnamed, then the first condition applies to the first variable, and the second condition applies to the second variable, and so on.
	MULTIPLE SETS OF CONDITIONS:
	A numeric matrix of conditioning values.
	Names are matched, same as in the simple use, except using column names.\cr
	Each row defines one set of conditions.

	Note that in univariate-conditional distributions, the number of conditioning variables needs to equal the total number of variables minus one.\cr
	By default (depending on the as.cset argument), a cks object is returned if there's one set of conditions, and a \link{cset} object is returned, if there's two or more sets conditions).}
\item{w}{Optional numeric vector of weights.}
\item{bw}{Postive numeric vector of length [1 or m], the bandwidth parameter(s).\cr
	If bw is missing, the bandwidth is computed for each variable using the bandwidth method (see bw.method below), which is subject to the smoothness parameter(s).}
\item{smoothness}{Positive numeric vector of length [1 or m], the relative smoothness parameter(s).\cr
	Ignored, if bw provided}
\item{kernel}{A (continuous) \link{kernel} object.}
\item{spline}{Logical, if true, use cubic Hermite splines as intermediate models.\cr In general, this should be true.}
\item{bw.method}{String, the bandwidth selection method.\cr Refer to \link{Bandwidth Selection}.}
\item{nc}{Integer, number of control points, in the spline.\cr Ignored, if spline is false.}
\item{tail}{
	Character vector of length [1 or M], either "lower" or "upper".\cr
	If lower (the default), lower tail probabilities, P (..., Xj <= xj, ...), are used.\cr
	If upper, upper tail probabilities, P (..., X >= xj, ...), are used.}
\item{trtype}{
	String, either "simple", "local" or "reflect", refer to vignette.\cr
	Note that the local method is used, where there's no trtype argument.}
\item{Xlim}{
	An ([1 or m] by 2) numeric matrix, giving the limits of each X variable.\cr
	The corresponding random variable is regarded as bounded, if either limit is finite.\cr
	In which case, a truncated smoothing algorithm is applied.}
\item{a, b}{
	Numeric vectors of length [1 or m].\cr
	This is an alternative way of specifying Xlim, above.}
\item{preserve.range}{Logical vector of length [1 or M].\cr
	If true, the default range used for range/sequence methods and plotting functions, will be the same as the original data.\cr
	If false (the default), the range is based on the data within a conditioning window, which is often smaller.}
\item{as.cset}{
	Logical, if true, a cset object is returned regardless of the number of conditions.\cr
	Note that this is ignored, if there is not one set of conditions.\cr\cr
	Also note that as.list.cset needs to be true, if as.cset is true, or the number of sets of conditions is not one.\cr
	(This requirement is temporary).}
\item{as.list.cset}{Logical, refer to as.cset, above.}
\item{warning}{Logical, if true, generate warning if there's no observations within the conditional window.}
\item{\dots}{Additional arguments not allowed}
}
\details{
PLEASE SET NOTES IN DESCRIPTION FIELD.

Here, M refers to the number of random variables (in the model), and m refers to the total number of variables.\cr
In nonconditional models, these are the same.\cr
In univariate-conditional models, M is one and m is equals to the number of conditions plus one.

Variable names use default values, if the x matrix is unnamed.\cr
(In conditional models, a warning is generated).

In conditional models, variables are reordered (internally), if there's named conditions.\cr
Note you can check names (and their order), with names method.\cr
(By default, only the conditional variables are returned).

Parameters that need to match m (the total number of variables), such as bw, smoothness and Xlim, will match the original data, regardless of the order of the conditions.\cr
i.e. The first bw value will match the first column in x, regardless of the order of the conditions.

However, parameters that need to match M (the number of random variables), such as tail, need to match the order of the resulting conditional variables.\cr
The same principle, applies to the x argument in the resulting function objects.

Refer to the vignette for more information.
}
\value{
Self-referencing function objects.

Refer to \link{Runtime Function Objects}

Except:\cr
The constructors for conditional distributions, return NULL, if there's no observations within the conditional window.\cr
(And by default, generate a warning).

They may also return cset objects, if as.cset is true, or there's more than one set of conditions.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Kernels}

\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}\cr
These can be used for statistical classification purposes.

\link{is.cks}, \link{ph.printf.phmodel}, \link{ph.plotf.cksuv}, \link{ph.plotf.cksmv}

\link{Bandwidth Selection}
}
\examples{
prep.ph.data ()

###########
#univariate
###########
uvfh <- pdfuv.cks (height)
uvFh <- cdfuv.cks (height)
uvFht <- qfuv.cks (height)

plot (uvfh)
plot (uvFh)
plot (uvfh, TRUE)

uvFht (c (0.25, 0.5, 0.75) )
ph.mode (uvfh)
ph.mode (uvfh, TRUE)

#############
#multivariate
#############
mvfh <- pdfmv.cks (trees2 [,-2])
mvFh <- cdfmv.cks (trees2 [,-2])

plot (mvfh, TRUE)
plot (mvfh,, TRUE)
plot (mvFh, TRUE)
plot (mvFh, fb = c (0.2, 0.8), clabs = c ("A", "B") )

##############
#conditional
#(single model)
###############
cons.c1 <- c (height=25, girth=35)
cfh <- pdfc.cks (trees2, conditions=cons.c1)
cFh <- cdfc.cks (trees2, conditions=cons.c1)
cFht <- qfc.cks (trees2, conditions=cons.c1)

plot (cfh)
plot (cFh)
plot (cfh, TRUE)

cFht (c (0.25, 0.5, 0.75) )
ph.mode (cfh)
ph.mode (cfh, TRUE)

#############
#conditional
#(two models)
#############
cons.c2 <- cbind (height=25:26, girth=35:36)

cfhs <- pdfc.cks (trees2, conditions=cons.c2, as.list.cset=TRUE)

plot (cfhs [[1]])

#########################
#multivariate-conditional
#########################
cons.mvc <- c (depth=311)

mvcFh <- cdfmvc.cks (quakes2, conditions=cons.mvc)
}
