\name{10_most_plot_methods}
\alias{plot.dkernel}
\alias{plot.ckernel}
\alias{plot.pmfuv.dks}
\alias{plot.cdfuv.dks}
\alias{plot.qfuv.dks}
\alias{plot.pdfuv.cks}
\alias{plot.cdfuv.cks}
\alias{plot.qfuv.cks}
\alias{plot.pdfc.cks}
\alias{plot.cdfc.cks}
\alias{plot.qfc.cks}
\alias{plot.chqf}
\alias{plot.pmfuv.cat}
\alias{plot.cdfuv.cat}
\alias{plot.qfuv.cat}
\alias{plot.cdf.el}
\alias{plot.qf.el}
\alias{plot.marginal.set}
\alias{plot.categorical.set}
\alias{plot.conditional.set}
\alias{lines.pdfuv}
\alias{lines.cdfuv}
\alias{lines.qfuv}
\title{Most Plots Methods}
\description{Plots of most objects in this package.}
\usage{
#kernels
\method{plot}{dkernel}(x, main, xlab="x", ylab,
    line.width=1, line.color="black", area.color, xlim, ylim, \dots,
	cdf=FALSE, bw=13)
\method{plot}{ckernel}(x, main, xlab="x", ylab,
    line.width=1, line.color="black", area.color, xlim = c (-1.05, 1.05), ylim, \dots,
    add=FALSE, cdf=FALSE)

#discrete kernel smoothing
\method{plot}{pmfuv.dks}(x, with.data.bars=FALSE, main, xlab, ylab="mass",
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{cdfuv.dks}(x, with.data.bars=FALSE, main, xlab, ylab="cumprob",
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{qfuv.dks}(x, no.data.bars, main, xlab="cumprob", ylab,
    line.width=1, line.color="black", area.color, \dots)

#continuous kernel smoothing
#(univariate)
\method{plot}{pdfuv.cks}(x, with.data.points=FALSE, main, xlab, ylab="density",
    line.width=1, line.color="black", area.color, \dots, add=FALSE)
\method{plot}{cdfuv.cks}(x, with.data.points=FALSE, main, xlab, ylab="cumprob",
    line.width=1, line.color="black", area.color, \dots, add=FALSE)
\method{plot}{qfuv.cks}(x, no.data.points, main, xlab="cumprob", ylab,
    line.width=1, line.color="black", area.color, \dots)

#continuous kernel smoothing
#(conditional)
\method{plot}{pdfc.cks}(x, no.data.points, main, xlab, ylab="density",
    line.width=1, line.color="black", area.color, \dots, add=FALSE)
\method{plot}{cdfc.cks}(x, no.data.points, main, xlab, ylab="cumprob",
    line.width=1, line.color="black", area.color, \dots, add=FALSE)
\method{plot}{qfc.cks}(x, no.data.points, main, xlab="cumprob", ylab,
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{chqf}(x, \dots)

#categorical distributions
\method{plot}{pmfuv.cat}(x, main, xlab, ylab="mass",
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{cdfuv.cat}(x, main, xlab, ylab="cumprob",
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{qfuv.cat}(x, main, xlab="cumprob", ylab,
    line.width=1, line.color="black", area.color, \dots)

#empirical-like distributions
\method{plot}{cdf.el}(x, with.data.points=FALSE, main, xlab, ylab="cumprob",
    line.width=1, line.color="black", area.color, \dots)
\method{plot}{qf.el}(x, with.data.points=FALSE, main, xlab="cumprob", ylab,
    line.width=1, line.color="black", area.color, \dots)

#distribution sets
\method{plot}{marginal.set}(x, main, xlab, ylab,
    colors, \dots)
\method{plot}{categorical.set}(x, main, xlab, ylab,
    legend=TRUE, colors, \dots)
\method{plot}{conditional.set}(x, main, xlab, ylab,
    legend=TRUE, colors, \dots)

\method{lines}{pdfuv}(x, \dots)
\method{lines}{cdfuv}(x, \dots)
\method{lines}{qfuv}(x, \dots)
}
\details{Refer to the vignette for more information.

Note that you can't plot chained quantile functions.
}
\arguments{
\item{x}{A suitable object.}
\item{with.data.bars}{If true, include a subpanel with data (frequency) bars.}
\item{with.data.points}{If true, include a subpanel with the data points.}
\item{no.data.bars}{Ignored.}
\item{no.data.points}{Ignored.}
\item{main}{.}
\item{xlab}{.}
\item{ylab}{.}
\item{xlim}{.}
\item{ylim}{.}
\item{line.width}{.}
\item{line.color}{.}
\item{area.color}{.}
\item{colors}{Colors for each plot.}
\item{add}{If true, add to an existing plot, otherwise, create a new plot.}
\item{bw}{Bandwidth.}
\item{cdf}{Plot the CDF rather than the PMF/PDF.}
\item{\dots}{Other arguments.}
\item{legend}{If true, add a legend to the plot.}
}
\examples{
data.prep ()

fh = pmfuv.dks (traffic.x, traffic.h, bw=25)
plot (fh, TRUE)
}
