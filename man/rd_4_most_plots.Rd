\name{rd_4_most_plots}
\alias{plot.nppdfuv}
\alias{plot.npcdfuv}
\alias{plot.npcdfuv.inverse}
\alias{plot.nppdfc}
\alias{plot.npcdfc}
\alias{plot.npcdfc.inverse}
\alias{plot.eccdf}
\alias{plot.eccdf.inverse}
\alias{plot.reglike}
\alias{lines.nppdfuv}
\alias{lines.npcdfuv}
\alias{lines.npcdfuv.inverse}
\alias{lines.nppdfc}
\alias{lines.npcdfc}
\alias{lines.npcdfc.inverse}
\alias{lines.eccdf}
\alias{lines.eccdf.inverse}
\alias{lines.reglike}
\title{Most Plots}
\description{Plots of most objects in this package.}
\usage{
#univariate
\method{plot}{nppdfuv}(x, with.points = FALSE, nmax = 2000, \dots)
\method{plot}{npcdfuv}(x, with.points = FALSE, nmax = 2000, \dots)
\method{plot}{npcdfuv.inverse}(x, \dots)
\method{lines}{nppdfuv}(x, \dots)
\method{lines}{npcdfuv}(x, \dots)
\method{lines}{npcdfuv.inverse}(x, \dots)

#conditional
\method{plot}{nppdfc}(x, \dots)
\method{plot}{npcdfc}(x, \dots)
\method{plot}{npcdfc.inverse}(x, \dots)
\method{lines}{nppdfc}(x, \dots)
\method{lines}{npcdfc}(x, \dots)
\method{lines}{npcdfc.inverse} (x, \dots)

#empirical
\method{plot}{eccdf}(x, with.points = FALSE, \dots)
\method{plot}{eccdf.inverse}(x, with.points = FALSE, \dots)
\method{lines}{eccdf}(x, \dots)
\method{lines}{eccdf.inverse}(x, \dots)

#regression-like
\method{plot}{reglike}(x, with.points = FALSE, \dots)
\method{lines}{reglike}(x, \dots)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{x}{An nppdfuv, npcdfuv, npcdfuv.inverse, nppdfc, npcdfc, npcdfc.inverse, eccdf, eccdf.inverse or reglike object.}
\item{with.points}{If true, plot the data points.}
\item{nmax}{If n <= nmax, plot all the data points. Otherwise, plot a random sample of the data points, with sample size = nmax.}
\item{\dots}{Other arguments.}
}
\examples{
data (trees)
nppdfuv.f = nppdfuv (trees$Height)
plot (nppdfuv.f, TRUE)
}
