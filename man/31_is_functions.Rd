\name{31_is_functions}
\alias{is.phob}
\alias{is.phmodel}
\alias{is.dpd}
\alias{is.cpd}
\alias{is.pmf}
\alias{is.dcdf}
\alias{is.dqf}
\alias{is.pdf}
\alias{is.ccdf}
\alias{is.cqf}
\alias{is.dpduv}
\alias{is.dpdc}
\alias{is.cpduv}
\alias{is.cpdmv}
\alias{is.cpdc}
\alias{is.cpdmvc}
\alias{is.pmfuv}
\alias{is.pmfc}
\alias{is.dcdfuv}
\alias{is.dcdfc}
\alias{is.dqfuv}
\alias{is.dqfc}
\alias{is.pdfuv}
\alias{is.pdfmv}
\alias{is.pdfc}
\alias{is.pdfmvc}
\alias{is.ccdfuv}
\alias{is.ccdfmv}
\alias{is.ccdfc}
\alias{is.ccdfmvc}
\alias{is.cqfuv}
\alias{is.cqfc}
\alias{is.cchqf}
\alias{is.dks}
\alias{is.cks}
\alias{is.cat}
\alias{is.el}
\alias{is.phspline}
\title{Is Functions}
\description{Functions to test if an object is of a particular class.}
\usage{
is.phob (xf)
is.phmodel (xf)

is.dpd (xf)
is.cpd (xf)

is.pmf (xf)
is.dcdf (xf)
is.dqf (xf)
is.pdf (xf)
is.ccdf (xf)
is.cqf (xf)

is.dpduv (xf, include.conditional=TRUE)
is.dpdc (xf, include.multivariate=TRUE)
is.cpduv (xf, include.conditional=TRUE)
is.cpdmv (xf, include.conditional=TRUE)
is.cpdc (xf, include.multivariate=TRUE)
is.cpdmvc (xf)

is.pmfuv (xf, include.conditional=TRUE)
is.pmfc (xf, include.multivariate=TRUE)
is.dcdfuv (xf, include.conditional=TRUE)
is.dcdfc (xf, include.multivariate=TRUE)
is.dqfuv (xf, include.conditional=TRUE)
is.dqfc (xf)

is.pdfuv (xf, include.conditional=TRUE)
is.pdfmv (xf, include.conditional=TRUE)
is.pdfc (xf, include.multivariate=TRUE)
is.pdfmvc (xf)
is.ccdfuv (xf, include.conditional=TRUE)
is.ccdfmv (xf, include.conditional=TRUE)
is.ccdfc (xf, include.multivariate=TRUE)
is.ccdfmvc (xf)
is.cqfuv (xf, include.conditional=TRUE)
is.cqfc (xf)

is.cchqf (xf)

is.dks (xf)
is.cks (xf)
is.cat (xf)
is.el (xf)

is.phspline (xf)
}
\arguments{
\item{xf}{An object to test.}
\item{include.conditional}{Logical, if true (the default), include conditional versions.}
\item{include.multivariate}{Logical, if true (the default), include multivariate versions.}
}
\details{
Note that DPD and CPD stand for discrete and continuous probability distributions, respectively.

A leading "d" is discrete and a leading "c" is continuous.
}
\value{A single logical value.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
ph.data.prep ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq, lower=0)

is.dpd (dfh)
is.cpd (dfh)
is.pmf (dfh)
is.dpduv (dfh)
}
