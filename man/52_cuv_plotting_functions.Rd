\name{52_cuv_plotting_functions}
\alias{plot_cpd}
\title{Plots of Continuous Univariate Models}
\description{Plots of continuous univariate probability distributions.}
\usage{
plot_cpd (sf, data=FALSE, ...,
    main, xlab, ylab,
    xlim, ylim,
    add=FALSE, axes=TRUE,
    line.width, line.color, fill.color)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this is a continuous univariate probability distribution.\cr\cr Refer to the references and see also sections.}
\item{data}{Logical, if true, include a subpanel with data points.\cr Ignored, if sf is a quantile function, or a conditional distribution.}
\item{main, xlab, ylab}{Optional strings, main/axes titles.}
\item{xlim, ylim}{Optional length-2 numeric vectors, giving the plot ranges.}
\item{add}{Logical, if true, add to an existing plot.}
\item{axes}{Logical vector of length one or two, if true, plot axis ticks with labels.}
\item{line.width}{Optional numeric, giving the main line width.\cr If missing, determined by global options.}
\item{line.color, fill.color}{Optional (R color) strings, giving the main line color and main fill color.\cr If missing, determined by global options.}
\item{\dots}{Ignored.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{set.ph.options}

\link{plot.cksuv}, \link{plot.eluv}

\link{plot_cpd_bv}, \link{plot_cpd_tv}
}
\examples{
ph.data.prep ()

cfh <- pdfuv.cks (height)

plot (cfh)
plot (cfh, TRUE)
}
