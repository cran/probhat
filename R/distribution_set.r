#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

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

marginal.set = function (constructor, x, ..., bw, smoothness)
{	m = ncol (x)
	pds = vector ("list", m)
	if (! missing (bw) )
	{	bw = .val.params (m, bw)
		for (j in 1:m)
			pds [[j]] = constructor (COL.of (x, j), ..., bw = bw [j])
	}
	else if (! missing (smoothness) )
	{	smoothness = .val.params (m, smoothness)
		for (j in 1:m)
			pds [[j]] = constructor (COL.of (x, j), ..., smoothness = smoothness [j])
	}
	else
	{	for (j in 1:m)
			pds [[j]] = constructor (COL.of (x, j), ...)
	}
	.distribution.set (m, pds, "marginal.set") 
}

categorical.set = function (constructor, x, ..., group.by)
{	group.by = as.character (group.by)
	ss = as.character (levels (as.factor (group.by) ) )
	m = length (ss)
	pds = vector ("list", m)
	for (j in 1:m)
	{	xsub = COL (x [group.by == ss [j] ], ss [j])
		pds [[j]] = constructor (xsub, ...)
	}
	.distribution.set.2 (m, pds, "categorical.set", .varname (x), ss) 
}

conditional.set = function (constructor, x, ..., group.by)
{	if (! is.matrix (group.by) )
		group.by = rbind (group.by)
	if (is.null (colnames (group.by) ) )
		labs = paste ("x", 1:ncol (group.by), sep="")
	else
		labs = colnames (group.by)
	m = nrow (group.by)
	pds = vector ("list", m)
	levnames = character (m)
	for (j in 1:m)
	{	levnames [j] = paste (labs, format (group.by [j,], digits=2), sep="=", collapse=", ")
		pds [[j]] = constructor (x, ..., conditions = group.by [j,])
	}
	.distribution.set.2 (m, pds, "conditional.set", .varname (x), levnames)  
}
