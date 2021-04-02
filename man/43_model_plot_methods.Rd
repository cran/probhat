\name{43_model_plot_methods}
\alias{ph.plotf.dksuv}
\alias{ph.plotf.cksuv}
\alias{ph.plotf.cksmv}
\alias{ph.plotf.cksc}
\alias{ph.plotf.cksmvc}
\alias{ph.plotf.catuv}
\alias{ph.plotf.catc}
\alias{ph.plotf.gmix}
\alias{ph.plotf.xmix}
\alias{ph.plotf.eluv}
\alias{ph.linesf.cpduv}
\title{Model Plot Methods}
\description{Plots methods for models, excluding distribution sets.}
\usage{
####################################
#discrete kernel smoothing models
#(call plot_dpd)
####################################
\method{ph.plotf}{dksuv}(sf, data=FALSE, \dots)

####################################
#continuous kernel smoothing models
#(call plot_cpd, plot_cpd_bv or plot_cpd_tv)
####################################
\method{ph.plotf}{cksuv}(sf, data=FALSE, \dots)
\method{ph.plotf}{cksmv}(sf, in3d=FALSE, data=FALSE, \dots)
\method{ph.plotf}{cksc}(sf, \dots)
\method{ph.plotf}{cksmvc}(sf, in3d=FALSE, data=FALSE, \dots)

####################################
#categorical models
#(call plot_dpd)
####################################
\method{ph.plotf}{catuv}(sf, \dots)
\method{ph.plotf}{catc}(sf, \dots)

####################################
#mixed input
####################################
\method{ph.plotf}{gmix}(sf, \dots)
\method{ph.plotf}{xmix}(sf, \dots)

####################################
#empirical-like models
#(call plot_cpd)
####################################
\method{ph.plotf}{eluv}(sf, data=FALSE, \dots)

####################################
#all continuous univariate models
####################################
\method{ph.linesf}{cpduv}(sf, \dots, xlim, n=200)
}
\arguments{
\item{sf}{A probability distribution.\cr Refer to the references and see also sections.}
\item{in3d}{Logical, if true, create a 3D plot.\cr Ignored, if sf has three or more random variables.}
\item{data}{If true, include a subpanel with data bars/points.\cr Ignored, if x is a quantile function, a conditional distribution, or has three or more random variables}
\item{xlim}{
	Length two numeric vector, giving plot range.\cr
	Currently, ignored for quantile functions.}
\item{n}{Integer, number of points.}
\item{\dots}{Other arguments for plot_dpd, plot_cpd, plot_cpd_bv and plot_cpd_tv.}
}
\details{
Refer to the vignette for more information.

Note that these methods call the functions plot_dpd, plot_cpd, plot_cpd_bv and plot_cpd_tv.\cr
Please refer to these functions for more information.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{plot_dpd}, \link{plot_cpd}

\link{plot_cpd_bv}, \link{plot_cpd_tv}
}
\examples{
prep.ph.data ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq)
cfh <- pdfuv.cks (height)
cfh2 <- pdfmv.cks (trees2 [,-2])

plot (dfh, TRUE)
plot (cfh, TRUE)
plot (cfh2,, TRUE)
}
