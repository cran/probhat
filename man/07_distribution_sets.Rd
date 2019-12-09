\name{07_distribution_sets}
\alias{marginal.set}
\alias{categorical.set}
\alias{conditional.set}
\title{Distribution Sets}
\description{Sets of distributions.}
\usage{
marginal.set (constructor, x, \dots, bw, smoothness)
categorical.set (constructor, x, \dots, group.by)
conditional.set (constructor, x, \dots, group.by)
}
\details{
These functions construct distribution set objects, which are a list of probability distributions.

Marginal sets:

The secondary constructor is applied to each column of x.


Categorical sets:

The secondary constructor is applied to subsets of (univariate) x, determined by group.by, which should be an integer-valued or character vector, same as x in pmfuv.cat.

Conditional sets:

The secondary constructor is applied to (multivariate) x. The group.by argument is a matrix of conditions, optionally with column names. Each row represents one set of conditions, with names taken from the column names, where applicable.

Note that currently, these functions (and their plots) may change significantly in future, and they have had relatively little error checking.
}
\arguments{
\item{constructor}{A univariate constructor (for marginal and categorical sets) or conditional constructor (for conditional sets).}
\item{x}{Vector (for categorical sets) or matrix (for marginal and conditional sets), corresponding to the x argument for the constructor used.}
\item{\dots}{Other arguments for the constructor, excluding conditions.}
\item{group.by}{Refer to details.}
\item{bw}{Bandwidth parameter.}
\item{smoothness}{Smoothness parameter.}
}
\examples{
data.prep ()

marset = marginal.set (qf.el, trees)
catset = categorical.set (pdfuv.cks, iris.Sepal.Length, group.by=iris.Species)

mconditions = cbind (Girth = marset [[1]](0.5),
    Height = marset [[2]](c (0.25, 0.5, 0.75)) )

conset = conditional.set (pdfc.cks, trees, group.by=mconditions)
}
