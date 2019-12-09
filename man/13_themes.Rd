\name{13_themes}
\alias{use.ph.theme}
\title{Themes}
\description{Set the color theme.}
\usage{
use.ph.theme (theme="blue", area.color, area.color.2="grey75", palette="Earth")
}
\arguments{
\item{theme}{Character, either "blue" or "green".}
\item{area.color}{Color for the area under univariate PMFs/PDFs/CDFs.}
\item{area.color.2}{Color for the area under univariate QFs, and data bars/points.}
\item{palette}{Color palette, for grDevices::hcl.colors, used when plotting distribution sets, and by the kernel.array function.}
}
\examples{
use.ph.theme ("blue")
}
