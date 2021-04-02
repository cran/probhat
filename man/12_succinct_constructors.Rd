\name{12_succinct_constructors}
\alias{Succinct Constructors}
\alias{pmf.dks}
\alias{cdf.dks}
\alias{qf.dks}
\alias{pdf.cks}
\alias{cdf.cks}
\alias{qf.cks}
\alias{pmf.cat}
\alias{cdf.cat}
\alias{qf.cat}
\alias{cdf.el}
\alias{qf.el}
\title{Succinct Constructors}
\description{Currently, these functions call the corresponding univariate constructors.}
\usage{
pmf.dks (\dots)
cdf.dks (\dots)
qf.dks (\dots)

pdf.cks (\dots)
cdf.cks (\dots)
qf.cks (\dots)

pmf.cat (\dots)
cdf.cat (\dots)
qf.cat (\dots)

cdf.el (\dots)
qf.el (\dots)
}
\arguments{
\item{\dots}{Argument list for the corresponding constructor.\cr Refer to the details section.}
}
\details{
Currently, these functions call the corresponding constructor with a uv suffix.

i.e.\cr
pmf.dks calls pmfuv.dks\cr
pdf.cks calls pdfuv.cks}
\value{Refer to univariate constructors.\cr i.e. pmfuv.dks.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}\cr
\link{Categorical Distributions}, \link{Empirical-Like Distributions}

\link{Conditional Distributions with Mixed Input Types}\cr
These can be used for statistical classification purposes.
}
\examples{
prep.ph.data ()

cFht <- qf.cks (height)
}
