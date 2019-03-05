\name{rd_6_plots_marginal}
\alias{plot.marginal}
\title{Plots of Marginal Distributions}
\description{Plots a marginal object.}
\usage{
\method{plot}{marginal}(x, nrow, ncol, with.points = FALSE, \dots)
}
\details{Refer to the vignette for more information. Note that nrow and ncol may be omitted. However, the plot layout may be suboptimal.}
\arguments{
\item{x}{A marginal object.}
\item{nrow}{Number of rows.}
\item{ncol}{Number of columns.}
\item{with.points}{If true, plot the data points.}
\item{\dots}{Other arguments.}
}
\examples{
data (trees)
trees = as.matrix (trees)
nppdfmv.f = nppdfmv (trees)
muv = marginal (nppdfmv.f)
plot (muv, 3, 1, TRUE)
}
