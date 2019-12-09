\name{11_bivariate_plot_methods}
\alias{plot.pdfmv.cks}
\alias{plot.cdfmv.cks}
\alias{plot.pdfmvc.cks}
\alias{plot.cdfmvc.cks}
\title{Plots of Bivariate Distributions}
\description{Plots of bivariate distributions.}
\usage{
#multivariate
\method{plot}{pdfmv.cks}(x, use.plot3d=FALSE,
    main, xlab, ylab, npoints=30, \dots, all=FALSE)
\method{plot}{cdfmv.cks}(x, use.plot3d=FALSE,
    main, xlab, ylab, npoints=30, \dots, all=FALSE)

#multivariate-conditional
\method{plot}{pdfmvc.cks}(x, use.plot3d=FALSE,
    main, xlab, ylab, npoints=30, \dots)
\method{plot}{cdfmvc.cks}(x, use.plot3d=FALSE,
    main, xlab, ylab, npoints=30, \dots)
}
\details{Refer to the vignette for more information.

Note that only multivariate probability distributions with two random variables can be plotted.

Also, note that these functions use the barsurf package.
}
\arguments{
\item{x}{A suitable object.}
\item{use.plot3d}{If true, use a 3D (surface) plot. Otherwise, use a 2D (contour/heatmap) plot.}
\item{main}{.}
\item{xlab}{.}
\item{ylab}{.}
\item{npoints}{The number of grid points in each direction.}
\item{all}{If true, plot both the contour plot and 3D plot, along with univariate plots and data.}
\item{\dots}{Other arguments.}
}
\examples{
data.prep ()

fh = pdfmv.cks (trees [,2:3])
plot (fh, TRUE)
}
