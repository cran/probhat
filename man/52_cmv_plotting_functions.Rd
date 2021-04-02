\name{52_cmv_plotting_functions}
\alias{plot_cpd_bv}
\alias{plot_cpd_tv}
\title{Plots of Continuous Multivariate Models}
\description{Plots of bivariate and trivariate continuous probability distributions.}
\usage{
#calls barsurf::plot_cfield or barsurf::plot_surface
plot_cpd_bv (sf, in3d=FALSE, data, \dots, n=30,
    main, xlab, ylab, zlab,
    xlim, ylim, zlim,
    add=FALSE, point.color)

#calls barsurf::plot_cfield_3d
plot_cpd_tv (sf, iso=FALSE, \dots,
    main, xlab, ylab, zlab,
    xlim, ylim, zlim,
    z.reverse=FALSE)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this is a continuous multivariate probability distribution, with two or three random variables.\cr\cr Refer to the references and see also sections.}
\item{in3d}{Logical, if true, produce a 3D plot.}
\item{iso}{Logical, if true, produce an isosurface plot.}
\item{data}{
	Logical, if true, plot data points.\cr
	Ignored, if in3d or add is true, or sf is a conditional distribution.\cr
	Defaults to true, if the number of observations is not more than 2000.}
\item{main, xlab, ylab, zlab}{Optional strings, main/axes titles.\cr\cr Note that these depend on the barsurf package.}
\item{xlim, ylim, zlim}{Optional length-2 numeric vectors, giving the plot ranges.}
\item{z.reverse}{Logical, if true, reverse the z axis.\cr Ignored, if zlim supplied.}
\item{n}{Numeric vector of length one or two, giving the number of grid points in each direction.}
\item{add}{
	Logical, if true, add to a previous plot.\cr
	Ignored, in 3d case.}
\item{point.color}{String, the (data) point color.}
\item{\dots}{Other arguments for barsurf functions, refer to details.}
}
\details{
These functions call barsurf::plot_cfield, barsurf::plot_surface, barsurf::plotf_cfield_3d and barsurf::plotf_isosurface.

There are some private intermediate functions, that change some of the defaults.\cr
Note that by default, more detail is used in bivariate CDFs than bivariate PDFs.
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{set.ph.options}

\link{plot_cpd}

\link{ph.plotf.cksmv}
}
\examples{
prep.ph.data ()

cfh2 <- pdfmv.cks (trees2 [,-2])

plot (cfh2)
plot (cfh2, TRUE)
plot (cfh2,, TRUE)
}
