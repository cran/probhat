\name{rd_2_nonparametric}
\alias{nppdfuv}
\alias{npcdfuv}
\alias{npcdfuv.inverse}
\alias{nppdfmv}
\alias{npcdfmv}
\alias{chained.npcdfmv.inverse}
\alias{nppdfc}
\alias{npcdfc}
\alias{npcdfc.inverse}
\alias{eccdf}
\alias{eccdf.inverse}
\alias{reglike}
\title{Nonparametric Models}
\description{Construct nonparametric probability distributions (that use kernel smoothing), empirical continuous cumulative distribution functions and their inverses, and regression-like models. Note that the objects are functions, which means that they can be evaluated.}
\usage{
#univariate
nppdfuv (x, spline = TRUE,
    kernel.pdf = sbc.pdf, kernel.cdf = NA,
    nc = 30, smoothness = 0.65, bw, w = NA)
npcdfuv (x, spline = TRUE,
    kernel.pdf = NA, kernel.cdf = sbc.cdf,
    nc = 30, smoothness = 0.65, bw, w = NA)
npcdfuv.inverse (x, kernel.pdf = NA, kernel.cdf = sbc.cdf,
    nc = 30, bw, smoothness = 0.65, w = NA)

#multivariate
nppdfmv (x, kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    smoothness = 0.65, bw, w = NA)
npcdfmv (x, kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    smoothness = 0.65, bw, w = NA)
chained.npcdfmv.inverse (
    x, kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    smoothness = 0.65, bw, w = NA)

#conditional
nppdfc (model, x, spline = TRUE,
    kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    nc = 30, smoothness = 0.65, bw, w = NA)
npcdfc (model, x, spline = TRUE,
    kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    nc = 30, smoothness = 0.65, bw, w = NA)
npcdfc.inverse (model, x, kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    nc = 30, smoothness = 0.65, bw, w = NA)

#empirical
eccdf (x, w = NA)
eccdf.inverse (x, w = NA)

#regression-like
reglike (model, x, kernel.pdf = sbc.pdf, kernel.cdf = sbc.cdf,
    nc = 30, nc.npc = nc, smoothness = 0.65, bw, w = NA)
}
\details{Refer to the vignette for more information.}
\arguments{
\item{model}{For conditional distributions, either a numeric vector giving the conditions, or a formula object. For regression-like models, a formula object. If it's a numeric vector, then it should have length m - 1, for the first m - 1 variables.}
\item{x}{For univariate distributions (including empirical continuous cumulative distribution functions and their inverses), x should be a numeric vector. Otherwise, x should be a matrix with n rows and m columns. Where rows represent data points (or observations) and columns represent variables.}
\item{spline}{If true (the default), use splines as intermediate models.}
\item{kernel.pdf}{A function giving the kernel PDF.}
\item{kernel.cdf}{A function giving the kernel CDF.}
\item{nc}{The number of control points in the spline.}
\item{nc.npc}{In regression-like models, the number of control points used in the conditional distributions.}
\item{smoothness}{A smoothness parameter of length one or m.}
\item{bw}{A bandwidth parameter of length one or m. If provided, the smoothness parameter is ignored. If missing, it's computed from the smoothness parameter.}
\item{w}{A numeric vector of weights that sum to approximately one.}
}
\examples{
data (trees)
nppdfuv.f = nppdfuv (trees$Height)
nppdfuv.f (80)
}
