#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.catuv = function (f, classes, g, h)
{	objs = .cat.data (g)
	gnames = nlevels = levels = n = m = g = NULL
	.UNPACK (objs)
 
	if (m != 1)
		stop ("uv model needs one variable")
	gname = gnames
	levels = levels [[1]]
	g = g [,1]

	h = .val.hvec (n, h)
	.gsum = sum (h)

	.probs = .iterate.uv (.pmfuv.cat.eval.scalar, n, g, h / .gsum, u=1:nlevels)
	.PROBS = cumsum (.probs)
	.PROBS [nlevels] = 1
	.EXTEND (f, classes,
		.probs, .PROBS, .gsum,
		gname, nlevels, glim = c (1, nlevels), levels, n, g, h)
}

.catc = function (f, classes, g, h, conditions, throw.warning)
{	objs = .cat.data (g, TRUE)
	gnames = nlevels = levels = n = m = g = NULL
	.UNPACK (objs)

	if (missing (conditions) )
		stop ("conditions required")
	ncon = length (conditions)
	if (ncon == 0)
		stop ("conditional models need at least one condition")
	M = m - ncon
	if (M != 1)
		stop ("uv conditional models need one random variable")
	h = .val.hvec (n, h)

	names = names (conditions)
	if (is.null (names) )
		names (conditions) = gnames [1:ncon]
	else
	{	J = match (names, gnames)
		if (any (is.na (J) ) )
			stop ("condition names not in variable names")
		J = c (J, (1:m)[-J])
		gnames = gnames [J]
		nlevels = nlevels [J]
		levels = levels [J]
		g = g [,J]
	}

	I = rep (TRUE, n)
	for (k in 1:ncon)
		I = I & .in.col (nlevels [k], levels [[k]], conditions [[k]], g [,k])
	
	n0 = n
	n = sum (I)
	if (n == 0)
	{	if (throw.warning)
			warning ("no observations within conditional window")
		NULL
	}
	else
	{	nlevels = nlevels [m]
		levels = levels [[m]]
		g = g [I, m, drop=FALSE]
		h = h [I]
		.gsum = sum (h)

		.probs = .iterate.uv (.pmfuv.cat.eval.scalar, n, g, h / .gsum, u=1:nlevels)
		.PROBS = cumsum (.probs)
		.PROBS [nlevels] = 1
		.EXTEND (f, classes,
			.probs, .PROBS, .gsum,
			gnames, conditions, nlevels, glim = c (1, nlevels), levels, n0, n, m, g, h)
	}
}

.in.col = function (nlevels, levels, condition, x)
{	condition = .cat.2int (nlevels, levels, condition, "unsuitable condition value")
	(condition == x)
}

.cat.2int = function (nlevels, levels, x, err)
{	if (is.integer (x) )
	{	if (any (x < 1 | x > nlevels) )
			stop (err)
		x
	}
	else if (is.character (x) )
	{	K = match (x, levels)
		if (any (is.na (K) ) )
			stop ("invalid category")
		.cat.2int (nlevels, levels, K, err)

	}
	else if (is.factor (x) )
		.cat.2int (nlevels, levels, as.character (x), err)
	else if (is.numeric (x) )
	{	y = as.integer (x)
		if (any (x != y) )
			stop (err)
		.cat.2int (nlevels, levels, y, err)
	}
	else
		stop (err)
}

#unpacked by .catuv, .catc, .mix
.cat.data = function (x, is.cond=FALSE)
{	if (is.matrix (x) )
	{	if (ncol (x) == 1)
		{	gname = colnames (x)
			if (is.null (gname) )
			{	warning ("applying default variable names, to single variable")
				gname = "g"
			}
			x = .cat.data.ext (x)
			x = .cat.data (list (g=x) )
			x$gnames = gname
			x
		}
		else
			stop ("currently, only 1-col matrices allowed")
	}
	else if (is.list (x) )
	{	m = length (x)
		if (m == 0)
			stop ("needs one or more categorical variables")
		gnames = .varnames.ext (m, names (x), "g", is.cond)

		y1 = .cat.data.ext (x [[1]])
		lev1 = levels (y1)

		n = length (y1)

		nlevels = integer (m)
		levels = vector ("list", m)
		y = matrix (0L, n, m)

		nlevels [1] = length (lev1)
		levels [[1]] = lev1
		y [,1] = as.integer (y1)

		if (m > 1)
		{	for (j in 2:m)
			{	yj = .cat.data.ext (x [[j]])
				levj = levels (yj)
				if (n != length (yj) )
					stop ("categorical input unequal lengths")

				nlevels [j] = length (levj)
				levels [[j]] = levj
				y [,j] = as.integer (yj)
			}
		}
		.LIST (gnames, nlevels, levels, n, m, g=y)
	}
	else
	{	x = .cat.data.ext (x)
		.cat.data (list (g=x) )
	}
}

.cat.data.ext = function (x)
{	if (is.integer (x) || is.character (x) )
	{	attributes (x) = NULL
		x = as.factor (x)
	}
	else if (is.factor (x) )
		NULL
	else if (is.numeric (x) )
	{	attributes (x) = NULL
		y = as.integer (x)
		if (any (x != y) )
			stop ("unsuitable numeric vector")
		x = as.factor (y)
	}
	else
		stop ("unsuitable object for categorical data")
	if (length (x) == 0)
		stop ("needs one or more values")
	if (any (is.na (x) ) )
		stop ("missing values")
	x
}

pmfuv.cat = function (g, h=1)
	.catuv (.pmfuv.cat.eval, .CV.pmfuv.cat, g, h)

cdfuv.cat = function (g, h=1)
	.catuv (.cdfuv.cat.eval, .CV.cdfuv.cat, g, h)

qfuv.cat = function (g, h=1)
	.catuv (.qfuv.cat.eval, .CV.qfuv.cat, g, h)

pmfc.cat = function (g, h=1, ..., conditions, warning=TRUE)
{	.arg.error (...)
	.catc (.pmfuv.cat.eval, .CV.pmfc.cat, g, h, conditions, warning)
}
	
cdfc.cat = function (g, h=1, ..., conditions, warning=TRUE)
{	.arg.error (...)
	.catc (.cdfuv.cat.eval, .CV.cdfc.cat, g, h, conditions, warning)
}

qfc.cat = function (g, h=1, ..., conditions, warning=TRUE)
{	.arg.error (...)
	.catc (.qfuv.cat.eval, .CV.qfc.cat, g, h, conditions, warning)
}

.pmfuv.cat.eval = function (g, ..., freq=FALSE, n)
{	. = .THAT ()
	x = .val.fg (.$nlevels, .$levels, g)
	p = .$.probs [x]
	.scale.freq (p, freq, .$.gsum, n)
}

.cdfuv.cat.eval = function (g, ..., freq=FALSE, n)
{	. = .THAT ()
	q = .val.fg (.$nlevels, .$levels, g)
	p = .$.PROBS [q]
	.scale.freq (p, freq, .$.gsum, n)
}

.qfuv.cat.eval = function (p, ..., level.names=FALSE)
{	. = .THAT ()
	.test.y.ok (p)
	x = .iterate.uv (.qfuv.cat.eval.2, .$.PROBS, u=p)
	if (level.names)
		.$levels [x]
	else
		x
}

.val.fg = function (nlevels, levels, x)
	.cat.2int (nlevels, levels, x, "unsuitable input for evaluation")

.pmfuv.cat.eval.scalar = function (n, x, h, u)
	sum (h [u == x])

.qfuv.cat.eval.2 = function (PROBS, y)
	which (y <= PROBS)[1]
