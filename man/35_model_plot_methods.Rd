\name{35_model_plot_methods}
\alias{plot.dksuv}
\alias{plot.cksuv}
\alias{plot.cksmv}
\alias{plot.cksc}
\alias{plot.cksmvc}
\alias{plot.catuv}
\alias{plot.eluv}
\alias{lines.cpduv}
\title{Model Plot Methods}
\description{Plots methods for models, excluding distribution sets.}
\usage{
####################################
#discrete kernel smoothing models
#(call plot_dpd)
####################################
\method{plot}{dksuv}(x, data=FALSE, \dots, freq)

####################################
#continuous kernel smoothing models
#(call plot_cpd, plot_cpd_bv or plot_cpd_tv)
####################################
\method{plot}{cksuv}(x, data=FALSE, \dots)
\method{plot}{cksmv}(x, in3d=FALSE, data=FALSE, \dots)
\method{plot}{cksc}(x, \dots)
\method{plot}{cksmvc}(x, in3d=FALSE, data=FALSE, \dots)

####################################
#categorical models
#(call plot_dpd)
####################################
\method{plot}{catuv}(x, \dots, combine, freq, space)

####################################
#empirical-like models
#(call plot_cpd)
####################################
\method{plot}{eluv}(x, data=FALSE, \dots)

####################################
#all continuous univariate models
####################################
\method{lines}{cpduv}(x, \dots)
}
\arguments{
\item{x}{A probability distribution.\cr Refer to the references and see also sections.}
\item{in3d}{Logical, if true, create a 3D plot.\cr Ignored, if sf has three or more random variables.}
\item{data}{If true, include a subpanel with data bars/points.\cr Ignored, if x is a quantile function, a conditional distribution, or has three or more random variables}
\item{combine}{Logical, if true, combine the bars.\cr By default, categorical distributions are separate and others are combined.}
\item{freq}{Logical, if true, return frequencies rather than probabilities.\cr Refer to the value section.\cr\cr Default depends on the object.\cr Refer to the corresponding constructors.}
\item{space}{Numeric, the space (in mm) between bars.\cr Ignored, unless combine is false.}
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
ph.data.prep ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq, lower=0)
cfh <- pdfuv.cks (height)

cfh2 <- pdfmv.cks (trees [,-2])
cfh3 <- pdfmv.cks (trees)

plot (dfh, TRUE)
plot (cfh, TRUE)
plot (cfh2,, TRUE)
plot (cfh3)
}
