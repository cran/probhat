\name{10_is_functions}
\alias{is.phob}
\alias{is.phpd}
\alias{is.phmodel}
\alias{is.dpd}
\alias{is.cpd}
\alias{is.pmf}
\alias{is.dcdf}
\alias{is.dqf}
\alias{is.pdf}
\alias{is.ccdf}
\alias{is.cqf}
\alias{is.pduv}
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
\alias{is.gmix}
\alias{is.xmix}
\alias{is.el}
\alias{is.phspline}
\title{Is Functions}
\description{Functions to test if an object is of a particular class.}
\usage{
is.phob (object)
is.phpd (object)
is.phmodel (object)

is.dpd (object)
is.cpd (object)

is.pmf (object)
is.dcdf (object)
is.dqf (object)
is.pdf (object)
is.ccdf (object)
is.cqf (object)

is.pduv (object, include.conditional=TRUE)

is.dpduv (object, include.conditional=TRUE)
is.dpdc (object, include.multivariate=TRUE)
is.cpduv (object, include.conditional=TRUE)
is.cpdmv (object, include.conditional=TRUE)
is.cpdc (object, include.multivariate=TRUE)
is.cpdmvc (object)

is.pmfuv (object, include.conditional=TRUE)
is.pmfc (object, include.multivariate=TRUE)
is.dcdfuv (object, include.conditional=TRUE)
is.dcdfc (object, include.multivariate=TRUE)
is.dqfuv (object, include.conditional=TRUE)
is.dqfc (object)

is.pdfuv (object, include.conditional=TRUE)
is.pdfmv (object, include.conditional=TRUE)
is.pdfc (object, include.multivariate=TRUE)
is.pdfmvc (object)
is.ccdfuv (object, include.conditional=TRUE)
is.ccdfmv (object, include.conditional=TRUE)
is.ccdfc (object, include.multivariate=TRUE)
is.ccdfmvc (object)
is.cqfuv (object, include.conditional=TRUE)
is.cqfc (object)

is.cchqf (object)

is.dks (object)
is.cks (object, include.xmix=TRUE)
is.cat (object, include.gmix=TRUE)
is.el (object)

is.gmix (object)
is.xmix (object)

is.phspline (object)
}
\arguments{
\item{object}{An object to test.}
\item{include.conditional}{Logical, if true (the default), include conditional versions.}
\item{include.multivariate}{Logical, if true (the default), include multivariate versions.}
\item{include.gmix}{Logical, if true (the default), include gmix objects.}
\item{include.xmix}{Logical, if true (the default), include xmix objects.}
}
\details{
Note that DPD and CPD stand for discrete and continuous probability distributions, respectively.

A leading "d" is discrete and a leading "c" is continuous.

Also, note that these relate to the objects, and not the number of variables.
(i.e. Object of class \link{pdfmv.cks}, are designed for multivariate models, but can be constructed with a single variable).
}
\value{A single logical value.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Succinct Constructors}\cr
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}
}
\examples{
prep.ph.data ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq)
is.dpd (dfh)
is.cpd (dfh)
}
