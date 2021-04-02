\name{90_global_options}
\alias{set.ph.options}
\alias{set.ph.theme}
\title{Global Options}
\description{Set the global options, including the color theme.}
\usage{
set.ph.options (\dots,
    rendering.style="r", theme="blue",
    main.line.width=1, main.line.color="#000000",
    main.fill.color, main.fill.color.2="#B0B0B0",
    semi.fill.color = "#00000030",
    palette="earth")

set.ph.theme (theme="blue")
}
\arguments{
\item{rendering.style}{Single character.\cr Refer to barsurf::set.bs.options.}
\item{theme}{Character, either "blue" or "green".\cr Other themes from the barsurf package are partially supported.}
\item{main.line.width}{Numeric, default main line width.}
\item{main.line.color}{String (R color string), default main line color.}
\item{main.fill.color}{String (R color string), default color for the area under univariate probability distributions, excluding quantile functions.}
\item{main.fill.color.2}{String (R color string), default color for the area under quantile functions, and univariate dks/cks data bars/points.}
\item{semi.fill.color}{String (R color string), default color for bivariate data points.}
\item{palette}{String, color palette, for grDevices::hcl.colors, used when plotting distribution sets, and by the kernel.array function.}
\item{\dots}{Ignored.}
}
\details{This function calls barsurf::set.bs.options, to set the rendering.style and theme.}
\references{Refer to the vignette for an overview, references and better examples.}
\examples{
set.ph.theme ("blue")
}
