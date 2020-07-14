#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.catuv = function (f, classes, g, w, freq=FALSE)
{	objs = .cat.data (g)
	variable.names = nlevels = levels = n = m = g = NULL
	UNPACK (objs)
 
	if (m != 1)
		stop ("uv model needs one variable")
	variable.name = variable.names
	levels = levels [[1]]
	g = g [,1]

	if (missing (w) ) w = NA
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w, FALSE)

	.probs = .iterate.uv (.pmfuv.cat.eval.scalar, is.weighted, n, g, w / sum (w), u=1:nlevels)
	.PROBS = cumsum (.probs)
	.PROBS [nlevels] = 1
	EXTEND (f, classes,
		.probs, .PROBS,
		is.weighted, freq,
		variable.name, nlevels, xlim = c (1, nlevels), levels, n, g, w)
}

.catc = function (f, classes, g, w, freq, conditions, throw.warning)
{	objs = .cat.data (g)
	variable.names = nlevels = levels = n = m = g = NULL
	UNPACK (objs)

	if (missing (conditions) )
		stop ("conditions required")
	ncon = length (conditions)
	if (ncon == 0)
		stop ("conditional models need at least one condition")
	M = m - ncon
	if (M != 1)
		stop ("uv conditional models need one random variable")

	names = names (conditions)
	if (is.null (names) )
		names (conditions) = variable.names [1:ncon]
	else
	{	J = match (names, variable.names)
		if (any (is.na (J) ) )
			stop ("condition names not in variable names")
		J = c (J, (1:m)[-J])
		variable.names = variable.names [J]
		nlevels = nlevels [J]
		levels = levels [J]
		g = g [,J]
	}
	if (missing (w) ) w = NA
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w, FALSE)
	
	I = rep (TRUE, n)
	for (k in 1:ncon)
		I = I & .in.col (nlevels [k], levels [[k]], conditions [[k]], g [,k])
	n = sum (I)
	if (n == 0)
	{	if (throw.warning)
			warning ("no observations within conditional window")
		NULL
	}
	else
	{	variable.name = variable.names [m]
		nlevels = nlevels [m]
		levels = levels [[m]]
		g = g [I, m, drop=FALSE]
		if (missing (w) ) w = NA
		if (is.weighted)
			w = w [I]
		.probs = .iterate.uv (.pmfuv.cat.eval.scalar, is.weighted, n, g, w / sum (w), u=1:nlevels)
		.PROBS = cumsum (.probs)
		.PROBS [nlevels] = 1
		EXTEND (f, classes,
			.probs, .PROBS,
			is.weighted, variable.name, conditions, freq, nlevels, xlim = c (1, nlevels), levels, n, g, w)
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
.cat.data = function (x)
{	if (is.list (x) )
	{	m = length (x)
		if (m == 0)
			stop ("needs one or more categorical variables")
		variable.names = .varnames.ext (m, names (x), "g")

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
		LIST (variable.names, nlevels, levels, n, m, g=y)
	}
	else
	{	x = .cat.data.ext (x)
		.cat.data (list (x=x) )
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

pmfuv.cat = function (g, h, ..., freq=FALSE)
	.catuv (.pmfuv.cat.eval, .CV.pmfuv.cat, g, h, freq)

cdfuv.cat = function (g, h, ...)
	.catuv (.cdfuv.cat.eval, .CV.cdfuv.cat, g, h)

qfuv.cat = function (g, h, ...)
	.catuv (.qfuv.cat.eval, .CV.qfuv.cat, g, h)

pmfc.cat = function (g, h, ..., conditions, warning=TRUE, freq=FALSE)
	.catc (.pmfuv.cat.eval, .CV.pmfc.cat, g, h, freq, conditions, warning)
	
cdfc.cat = function (g, h, ..., conditions, warning=TRUE)
	.catc (.cdfuv.cat.eval, .CV.cdfc.cat, g, h, FALSE, conditions, warning)

qfc.cat = function (g, h, ..., conditions, warning=TRUE)
	.catc (.qfuv.cat.eval, .CV.qfc.cat, g, h, FALSE, conditions, warning)

.pmfuv.cat.eval = function (g, ..., freq)
{	. = THAT ()
	x = .val.g (.$nlevels, .$levels, g)
	if (missing (freq) )
		freq = .$freq
	p = .$.probs [x]
	if (freq)
	{	if (.$is.weighted) 
			sum (.$w) * p
		else
			sum (.$n) * p
	}
	else
		p
}

.cdfuv.cat.eval = function (q)
{	. = THAT ()
	q = .val.g (.$nlevels, .$levels, q)
	.$.PROBS [q]
}

.qfuv.cat.eval = function (p, ..., name=FALSE)
{	. = THAT ()
	.test.y.ok (p)
	x = .iterate.uv (.qfuv.cat.eval.2, .$.PROBS, u=p)
	if (name)
		.$levels [x]
	else
		x
}

.val.g = function (nlevels, levels, x)
	.cat.2int (nlevels, levels, x, "unsuitable input for evaluation")

.pmfuv.cat.eval.scalar = function (is.weighted, n, x, w, u)
{	if (is.weighted)
		sum (w [u == x])
	else
		sum (u == x) / n
}

.qfuv.cat.eval.2 = function (PROBS, y)
	which (y <= PROBS)[1]
