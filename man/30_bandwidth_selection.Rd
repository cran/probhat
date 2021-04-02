\name{30_bandwidth_selection}
\alias{Bandwidth Selection}
\alias{auto.dbw}
\alias{auto.cbw}
\title{Bandwidth Selection}
\description{Functions for bandwidth selection.

THESE FUNCTIONS SHOULD BE REGARDED AS SUB-OPTIMAL.}
\usage{
auto.dbw (x, ..., bw.method="ph.default", smoothness=1)
auto.cbw (x, ..., bw.method="ph.default", smoothness=1)
}
\arguments{
\item{x}{A numeric vector, of data.\cr\cr In the continuous case, x may also be a matrix.}
\item{bw.method}{String, the initial bandwidth selection method.\cr Currently, "ph.default", "Scott" or "Silverman".\cr Refer to details.}
\item{smoothness}{Numeric, smoothness (scaling) parameter.\cr Refer to details.}
\item{\dots}{Ignored.}
}
\details{THESE FUNCTIONS SHOULD BE REGARDED AS SUB-OPTIMAL.

These functions computes an initial bandwidth.\cr
Then the initial bandwidth parameter is multiplied by the smoothness parameter.

In the discrete case (auto.dbw), this bandwidth is rounded up to the nearest odd integer.

Currently, there are three options:

(1) bw.method="ph.default".

For a single variable/column, the bandwidth is equal to the difference between the quantiles, marking the middle 0.66 of observations.\cr
i.e. diff (quantile (x, c (0.17, 0.83) ) )

For m variables/columns, 0.66 is replaced with 0.66^(1 / m).

(2) bw.method="Scott", which calls stats::bw.nrd, for each variable/column.

(3) bw.method="Silverman", which call stats::bw.nrd0, for each variable/column.
}
\value{In the discrete case, a single integer.\cr In the continuous case, a numeric vector.}
\references{
Refer to the vignette for an overview, references and better examples.

Also please refer to stats::bw.nrd and stats::bw.nrd0 for references, and more information.
}
\seealso{
\link{bw.nrd}, \link{bw.nrd0}

\link{Discrete Kernel Smoothing}, \link{Continuous Kernel Smoothing}
}
\examples{
prep.ph.data ()
auto.cbw (trees2)
}
