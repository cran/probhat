#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.distribution.set = function (npds, pds, class)
	EXTEND (pds, c (class, "distribution.set"), npds)

.distribution.set.2 = function (npds, pds, class, varname, levnames)
	EXTEND (pds, c (class, "distribution.set"), npds, varname, levnames)

.gset = function (constructor, group.by, x, ...)
{	if (is.list (group.by) )
		group.by = group.by [[1]]
	group.by = as.character (group.by)
	ss = as.character (levels (as.factor (group.by) ) )
	m = length (ss)
	pds = vector ("list", m)
	for (j in 1:m)
	{	xsub = x [group.by == ss [j] ]
		xsub = cbind (xsub)
		colnames (xsub) = ss [j]
		pds [[j]] = constructor (xsub, ...)
	}
	.distribution.set.2 (m, pds, "ph3.gset", .varname (x), ss) 
}

.mset.cks = function (constructor, x, bw, smoothness, ...)
{	m = ncol (x)
	pds = vector ("list", m)
	if (! missing (bw) )
	{	bw = .val.params (m, bw)
		for (j in 1:m)
			pds [[j]] = constructor (x [,j, drop=FALSE], ..., bw = bw [j])
	}
	else
	{	smoothness = .val.params (m, smoothness)
		for (j in 1:m)
			pds [[j]] = constructor (x [,j, drop=FALSE], ..., smoothness = smoothness [j])
	}
	.distribution.set (m, pds, "ph3.mset") 
}

.mset.el = function (constructor, x, ...)
{	m = ncol (x)
	pds = vector ("list", m)
	for (j in 1:m)
		pds [[j]] = constructor (x [,j, drop=FALSE], ...)
	.distribution.set (m, pds, "ph3.mset") 
}

pdfuv.gset.cks = function (x, ..., group.by) .gset (pdfuv.cks, group.by, x, ...)
cdfuv.gset.cks = function (x, ..., group.by) .gset (cdfuv.cks, group.by, x, ...)
qfuv.gset.cks = function (x, ..., group.by) .gset (qfuv.cks, group.by, x, ...)
cdfuv.gset.el = function (x, ..., group.by) .gset (cdfuv.el, group.by, x, ...)
qfuv.gset.el = function (x, ..., group.by) .gset (qfuv.el, group.by, x, ...)

pdfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (pdfuv.cks, x, bw, smoothness, ...)
cdfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (cdfuv.cks, x, bw, smoothness, ...)
qfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (qfuv.cks, x, bw, smoothness, ...)
cdfuv.mset.el = function (x, ...) .mset.el (cdfuv.el, x, ...)
qfuv.mset.el = function (x, ...) .mset.el (qfuv.el, x, ...)
