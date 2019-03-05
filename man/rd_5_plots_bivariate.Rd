\name{rd_5_plots_bivariate}
\alias{plot.nppdfmv}
\alias{plot.npcdfmv}
\alias{plot.chained.npcdfmv.inverse}
\title{Plots of Bivariate Distributions}
\description{Plots of bivariate nonparametric distributions.}
\usage{
\method{plot}{nppdfmv}(x, use.plot3d = FALSE,
    xlab, ylab, npoints = 30, \dots, all = FALSE)
\method{plot}{npcdfmv}(x, use.plot3d = FALSE,
    xlab, ylab, npoints = 30, \dots)
\method{plot}{chained.npcdfmv.inverse}(x, \dots)
}
\details{Refer to the vignette for more information. Note that only bivariate nppdfmv and bivariate npcdfmv objects (with m = 2) can be plotted. Other objects produce error messages.}
\arguments{
\item{x}{An nppdfmv or npcdfmv object.}
\item{use.plot3d}{If true, use a 3d plot. Otherwise, use a contour plot.}
\item{xlab}{The x axis label.}
\item{ylab}{The y axis label.}
\item{npoints}{The number of grid points in each direction.}
\item{all}{If true, plot all combinations.}
\item{\dots}{Other arguments.}
}
\examples{
data (trees)

#matrix with height and volume
trees.sub = as.matrix (trees)[,2:3]

nppdfmv.f = nppdfmv (trees.sub)
plot (nppdfmv.f)
}
