#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.distribution.set = function (npds, pds, class)
	.EXTEND (pds, c (class, "dset", "phob"), npds)

.distribution.set.2 = function (npds, pds, class, varname, levnames)
	.EXTEND (pds, c (class, "dset", "phob"), npds, varname, levnames)

.gset = function (constructor, group.by, x, ...)
{	if (is.list (group.by) )
		group.by = group.by [[1]]
	if (is.factor (group.by) )
	{	ss = as.character (levels (group.by) )
		group.by = as.character (group.by)
	}
	else
	{	group.by = as.character (group.by)
		ss = as.character (levels (as.factor (group.by) ) )
	}
	m = length (ss)
	pds = vector ("list", m)
	x = cbind (x)
	for (j in 1:m)
	{	xsub = x [group.by == ss [j],, drop=FALSE]
		pds [[j]] = constructor (xsub, ...)
	}
	.distribution.set.2 (m, pds, "ph4.gset", .varname (x), ss)
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
	.distribution.set (m, pds, "ph4.mset")
}

.mset.el = function (constructor, x, ...)
{	m = ncol (x)
	pds = vector ("list", m)
	for (j in 1:m)
		pds [[j]] = constructor (x [,j, drop=FALSE], ...)
	.distribution.set (m, pds, "ph4.mset") 
}

ph4.pdfuv.gset.cks = function (g, x, ...)
	.EXTEND (.gset (pdfuv.cks, g, x, ...), "ph4.pdfuv.gset.cks")
ph4.cdfuv.gset.cks = function (g, x, ...) .gset (cdfuv.cks, g, x, ...)
ph4.qfuv.gset.cks = function (g, x, ...) .gset (qfuv.cks, g, x, ...)
ph4.cdfuv.gset.el = function (g, x, ...) .gset (cdfuv.el, g, x, ...)
ph4.qfuv.gset.el = function (g, x, ...) .gset (qfuv.el, g, x, ...)

ph4.pdfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (pdfuv.cks, x, bw, smoothness, ...)
ph4.cdfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (cdfuv.cks, x, bw, smoothness, ...)
ph4.qfuv.mset.cks = function (x, ..., bw, smoothness=1) .mset.cks (qfuv.cks, x, bw, smoothness, ...)
ph4.cdfuv.mset.el = function (x, ...) .mset.el (cdfuv.el, x, ...)
ph4.qfuv.mset.el = function (x, ...) .mset.el (qfuv.el, x, ...)

ph4.pdfmv.gset.cks = function (g, x, ...)
	.EXTEND (.gset (pdfmv.cks, g, x, ...), "ph4.pdfmv.gset.cks")

as.list.dset = function (x, ...)
{	attributes (x) = NULL
	x
}
