\name{18_other_functions}
\alias{COL}
\alias{COL.of}
\alias{plot_2x2}
\alias{data.prep}
\title{Other Functions}
\description{Other functions, not covered in other help files.}
\usage{
COL (x, name="x")
COL.of (x, \dots)

plot_2x2 (f1, f2, f3, f4, main.1, main.2, main.3, main.4, \dots)

data.prep (eval=TRUE, echo=FALSE)
}
\details{
The COL and COL.of functions are designed to create a single-column matrix, with a column name.

The COL.of function is equivalent to column subsetting with drop=FALSE.

The plot_2x2 function, plots a two-by-two array of objects, closely spaced, without xlab and ylab labels.\cr
Note that this function may be removed, in the future.

The data.prep function, which is only designed for tests and examples, runs an R script stored inside the function, itself.
}
\arguments{
\item{x}{For COL, a vector. For COL.of, a matrix, preferably with column names, and optionally with row names.}
\item{name}{Name of the variable.}
\item{\dots}{For COL.of, the column indices.}
\item{f1, f2, f3, f4}{Objects.}
\item{main.1, main.2, main.3, main.4}{Main titles.}
\item{eval}{Run the script.}
\item{echo}{Echo the script.}
}
\examples{
data.prep ()

COL.of (trees, 2)
}
