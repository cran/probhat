\name{53_cmv_plotting_functions}
\alias{plot_cpd_bv}
\alias{plot_cpd_tv}
\title{Plots of Continuous Multivariate Models}
\description{Plots of bivariate and trivariate continuous probability distributions.}
\usage{
#calls barsurf::plot_cfield or barsurf::plot_surface
plot_cpd_bv (sf, in3d=FALSE, data=FALSE, \dots, all=FALSE, n=30,
    main, xlab, ylab,
    xlim, ylim, zlim)

#calls barsurf::plot_cfield_3d
plot_cpd_tv (sf, \dots,
    main, xlab, ylab, zlab,
    xlim, ylim, zlim,
    reverse.z=FALSE)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this is a continuous multivariate probability distribution, with two or three random variables.\cr\cr Refer to the references and see also sections.}
\item{in3d}{Logical, if true, produce a 3D plot.}
\item{data}{Logical, if true, plot data points.\cr Ignored, if in3d is true, or sf is a conditional distribution.}
\item{main, xlab, ylab, zlab}{Optional strings, main/axes titles.\cr\cr Note that these depend on the barsurf package.\cr And at the time of writing (with barsurf version 0.5.0), zlab is not supported.}
\item{xlim, ylim, zlim}{Optional length-2 numeric vectors, giving the plot ranges.}
\item{reverse.z}{Logical, if true, reverse the z axis.\cr Ignored, if zlim supplied.}
\item{all}{Logical, if true, plot a 2x2 plot array, with different subplots.\cr Ignored, if sf is a conditional distribution.}
\item{n}{Numeric vector of length one or two, giving the number of grid points in each direction.}
\item{\dots}{Other arguments for barsurf::plot_cfield, barsurf::plot_surface, plot::plot_cfield_3d.}
}
\details{These functions call barsurf::plot_cfield, barsurf::plot_surface, plot::plot_cfield_3d.}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{set.ph.options}

\link{plot_cpd}

\link{plot.cksmv}
}
\examples{
ph.data.prep ()

cfh2 <- pdfmv.cks (trees [,-2])
cfh3 <- pdfmv.cks (trees)

plot (cfh2)
plot (cfh2, TRUE)
plot (cfh2,, TRUE)
plot (cfh2, all=TRUE)

plot (cfh3)
}
