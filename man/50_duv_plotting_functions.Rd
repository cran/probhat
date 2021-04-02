\name{50_duv_plotting_functions}
\alias{plot_dpd}
\title{Plots of Discrete Univariate Models}
\description{Plots of discrete univariate probability distributions.}
\usage{
plot_dpd (sf, data=FALSE, ...,
    main, xlab, ylab,
    xlim, ylim,
    add=FALSE, axes=TRUE,
    combine = is.dks (sf), freq=FALSE, n, space=0,
    line.width, line.color, fill.color)
}
\arguments{
\item{sf}{A suitable function object.\cr Here, this is a discrete univariate probability distribution.\cr\cr Refer to the references and see also sections.}
\item{data}{Logical, if true, include a subpanel with the data bars.\cr Ignored, if sf is a quantile function or a conditional distribution.}
\item{main, xlab, ylab}{Optional strings, main/axes titles.}
\item{xlim, ylim}{Optional length-2 numeric vectors, giving the plot ranges.}
\item{add}{Logical, if true, add to an existing plot.}
\item{axes}{Logical vector of length one or two, if true, plot axis ticks with labels.}
\item{combine}{Logical, if true, combine the bars.}
\item{freq}{
	Logical, if true, plot frequencies.\cr
	Currently, ignored in quantile functions.}
\item{n}{
	An integer, representing the sample size, and is ignored except for PMFs/CDFs with freq set to TRUE.\cr
	(If missing, it defaults to the number of observations used, or the sum of their frequencies).}
\item{space}{Numeric, the space (in mm) between the bars.\cr Ignored, if combine is true.}
\item{line.width}{Optional numeric, giving the main line width.\cr If missing, determined by global options.}
\item{line.color, fill.color}{Optional (R color) strings, giving the main line color and main fill color.\cr If missing, determined by global options.}
\item{\dots}{Ignored.}
}
\references{Refer to the vignette for an overview, references and better examples.}
\seealso{
\link{set.ph.options}

\link{ph.plotf.dksuv}, \link{ph.plotf.catuv}
}
\examples{
prep.ph.data ()

dfh <- pmfuv.dks (traffic.bins, traffic.freq)

plot (dfh)
plot (dfh, TRUE)
}
