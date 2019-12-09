\name{08_is_functions}
\alias{is.dpd}
\alias{is.cpd}
\alias{is.pmfuv}
\alias{is.pdfuv}
\alias{is.cdfuv}
\alias{is.qfuv}
\alias{is.chqf}
\title{Is Functions}
\description{
Functions to test if an object is of a particular (super) class.

Note that DPD and CPD stand for discrete and continuous probability distributions, respectively.
}
\usage{
is.dpd (f)
is.cpd (f)

is.pmfuv (f)
is.pdfuv (f)
is.cdfuv (f)
is.qfuv (f)

is.chqf (f)
}
\arguments{
\item{f}{The object to test.}
}
\examples{
data.prep ()

fh = pmfuv.dks (traffic.x, traffic.h, bw=25)
is.pmfuv (fh)
}
