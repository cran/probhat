#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.n.unique = function (x) length (unique (x) )
.any.duplicates = function (x) .n.unique (x) != length (x)

auto.dbw = function (x, ..., bw.method="ph.default", smoothness=1)
{	bw = auto.cbw (x, ..., bw.method=bw.method, smoothness=smoothness)
	bw = as.integer (round (bw) )
	if (bw %% 2 == 0)
		bw = bw + 1
	bw
}

auto.cbw = function (x, ..., bw.method="ph.default", smoothness=1)
{	bw.method = tolower (bw.method)
	if (! is.matrix (x) )
		x = cbind (x)
	if (bw.method == "ph.default")
	{	m = ncol (x)
		P = 0.66^(1 / m)
		a = (1 - P) / 2
		b = 1 - a
		bw = numeric (m)
		for (j in seq_len (m) )
			bw [j] = diff (quantile (x [,j], c (a, b) ) )
	}
	else if (bw.method == "scott") bw = apply (x, 2, bw.nrd)
	else if (bw.method == "silverman") bw = apply (x, 2, bw.nrd0)
	else
		stop ("bw.method needs to be ph.default, Scott or Silverman")
	smoothness * bw
}

.midpoints = function (x)
{	n = length (x)
	(x [-n] + x [-1]) / 2
}

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

.varname = function (x)
{	if (is.matrix (x) && ncol (x) == 1)
		colnames (x)
	else
		"x"
}

.varnames = function (x)
{	if (is.matrix (x) )
		.varnames.ext (ncol (x), colnames (x) )
	else
		"x"
}

.varnames.ext = function (m, variable.names, prefix="x")
{	if (is.null (variable.names) )
		paste (prefix, 1:m, sep="")
	else
	{	if (m != .n.unique (variable.names) )
			stop ("needs unique variable names")
		if (any (variable.names == "") )
			stop ("needs all (or no) variable names")
		variable.names
	}
}

.blabs = function (x)
{	if (is.matrix (x) )
		rownames (x)
	else
		names (x)
}

.val.k = function (pkernel, ..., super="ckernel")
{	if (is.function (pkernel) ) kernel = pkernel (...)
	else kernel = 0

	if (inherits (kernel, super) )
		kernel
	else
		stop ("needs constructor for suitable kernel object")
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

.val.x.uv = function (x, one.or.more=FALSE)
{	attributes (x) = NULL
	if (is.vector (x) )
	{	x = as.numeric (x)
		if (one.or.more && length (x) == 0)
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

.val.x.uv.or.mv = function (x)
{	if (is.matrix (x) )
		.val.x.mv (x)
	else
		cbind (.val.x.uv (x) )
}

.val.w = function (is.weighted, n, w, scale=TRUE)
{	if (is.weighted)
	{	w = as.numeric (w)
		if (n != length (w) )
			stop ("length (w) != number of observations")
		if (any (! is.finite (w) ) )
			stop ("all w values need to be finite")
		if (any (w <= 0) )
			stop ("all w value need to be >= 0")
		if (scale)
			w = w / sum (w)
		w
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

.deflab = function (f, lab)
{	if (missing (lab) )
	{	if (is.dpdc (f) || is.cpdc (f) )
			paste ("conditional", f %$% "variable.name")
		else
			f %$% "variable.name"
	}
	else
		lab
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
