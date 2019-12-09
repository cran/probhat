#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.as.integer.matrix = function (x)
{	y = as.integer (x)
	dim (y) = dim (x)
	y
}

.as.numeric.matrix = function (x)
{	y = as.numeric (x)
	dim (y) = dim (x)
	y
}

COL = function (x, name="x")
{	x = matrix (x,, 1)
	colnames (x) = name
	x
}

COL.of = function (x, ...)
	x [, c (...), drop=FALSE]

.varname = function (x)
{	if (is.matrix (x) && ncol (x) == 1)
		colnames (x)
	else
		"x"
}

.varnames = function (x)
{	m = ncol (x)
	variable.names = colnames (x)
	if (is.null (variable.names) )
		paste ("x", 1:m, sep="")
	else
	{	if (m != .n.unique (variable.names) )
			stop ("x needs unique variable names")
		variable.names
	}
}

.blabs = function (x)
{	if (is.matrix (x) )
		rownames (x)
	else
		names (x)
}

.val.k = function (kernel, ...)
{	if (is.function (kernel) )
		kernel (...)
	else if (inherits (kernel, "kernel") )
		kernel
	else
		stop ("needs kernel object, or suitable constructor")
}

.val.params = function (m, param)
{	nparam = length (param)
	if (nparam == 1)
		rep (param, m)
	else if (nparam == m)
		param
	else
		stop ("parameter needs to have length 1 or m")
}

.val.x.uv = function (x)
{	attributes (x) = NULL
	if (is.vector (x) )
	{	x = as.numeric (x)
		if (length (x) < 1)
			stop ("x needs one or more values")
		if (any (! is.finite (x) ) )
			stop ("all x values need to be finite")
		x
	}
	else
		stop ("needs vector (or matrix)")
}

.val.x.mv = function (x)
{	if (! is.matrix (x) )
		stop ("multivariate models need matrix")
	x = .as.numeric.matrix (x)
	if (nrow (x) < 1)
		stop ("x needs one or more rows")
	if (any (! is.finite (x) ) )
		stop ("all x values need to be finite")
	x
}

.val.w = function (is.weighted, n, w)
{	if (is.weighted)
	{	w = as.numeric (w)
		if (n != length (w) )
			stop ("length (w) != number of observations")
		if (any (! is.finite (w) ) )
			stop ("all w values need to be finite")
		if (any (w <= 0) )
			stop ("all w value need to be >= 0")
		w / sum (w)
	}
	else
		NA
	
}

.val.u.mv = function (m, u)
{	if (! is.matrix (u) )
		u = rbind (u)
	if (m != ncol (u) )
		stop ("incorrect number of columns")
	u
}

.iterate.uv = function (f, ..., u)
{	n = length (u)
	y = numeric (n)
	for (i in seq_len (n) )
		y [i] = f (..., u [i])
	y
}

.iterate.mv = function (f, ..., u)
{	n = nrow (u)
	y = numeric (n)
	for (i in seq_len (n) )
		y [i] = f (..., u [i,])
	y
}

.iterate.mv.2 = function (f, ..., y)
{	n = nrow (y)
	x = numeric (n)
	for (i in seq_len (n) )
		x [i,] = f (..., y [i,])
	x
}

.test.y.ok = function (y)
{	if (any (y < 0 | y > 1) )
		stop ("probabilities need to be between 0 and 1")
}
