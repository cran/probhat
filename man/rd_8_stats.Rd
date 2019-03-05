\name{rd_8_stats}
\alias{npmean}
\alias{npmode}
\alias{nprng}
\title{Statistics and Random Number Generation}
\description{Compute the mean and mode from nonparametric PDFs, including conditional distributions. Also, generate random numbers (or random samples) using inversion sampling.}
\usage{
npmean (nppdfuv.f)
npmode (nppdfuv.f, include.boundaries = TRUE,
    all = FALSE, warning = FALSE)
nprng (npcdf.f.inverse, n)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{nppdfuv.f}{An nppdfuv object.}
\item{npcdf.f.inverse}{An npcdfuv.inverse or chained.npcdfmv.inverse object.}
\item{include.boundaries}{Include the boundary points, if they're modal points, which isn't possible without modifying the splines.}
\item{all}{Return all modal points.}
\item{warning}{Produce a warning, if zero or multiple modal points.}
\item{n}{Number of random numbers.}
}
\examples{
data (trees)
npmode (nppdfuv (trees$Height) )
}
