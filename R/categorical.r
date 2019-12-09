#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.int.cat = function (x)
{	if (is.matrix (x) )
		y = .as.integer.matrix (x)
	else
		y = as.integer (x)
	if (! all (x == y) )
		.cat.err ()
	y
}

.cat.err = function ()
	stop ("x needs to be integer-valued or character")

.catuv = function (f, superclass, subclass, x, w, sample.space)
{	variable.name = .varname (x)
	attributes (x) = NULL
	n = length (x)
	if (is.integer (x) )
	{	if (missing (sample.space) )
		{	nlevels = .n.unique (x)
			sample.space = as.character (1:nlevels)
		}
		else
		{	sample.space = as.character (sample.space)
			nlevels = length (sample.space)
			if (nlevels != .n.unique (x) )
				stop ("incorrect number of levels")
		}
	}
	else if (is.character (x) )
	{	if (! missing (sample.space) )
			warning ("if x character, sample.space ignored")
		y = as.factor (x)
		sample.space = levels (y)
		nlevels = length (sample.space)
		x = as.integer (y)
	}
	else if (is.numeric (x) )
		.catuv (f, superclass, subclass, .int.cat (x), w, sample.space)
	else
		.cat.err ()
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w)
	.probs = .iterate.uv (.pmfuv.cat.eval.scalar, is.weighted, n, x, w, u=1:nlevels)
	.PROBS = cumsum (.probs)
	.PROBS [nlevels] = 1
	EXTEND (f, c (subclass, "cat", superclass, "dpd", "phmodel"),
		.probs, .PROBS,
		variable.name, is.weighted, nlevels, levels=sample.space, n, x, w)
}

pmfuv.cat = function (x, sequence, w=NA, levels)
{	f = function (x)
	{	. = THAT ()
		x = .val.cat.x.uv (.$nlevels, .$levels, x)
		.$.probs [x]
	}
	.catuv (f, "pmfuv", "pmfuv.cat", x, w, levels)
}

cdfuv.cat = function (x, sequence, w=NA, levels)
{	f = function (x)
	{	. = THAT ()
		x = .val.cat.x.uv (.$nlevels, .$levels, x)
		.$.PROBS [x]
	}
	.catuv (f, "cdfuv", "cdfuv.cat", x, w, levels)
}

qfuv.cat = function (x, sequence, w=NA, levels)
{	f = function (y)
	{	. = THAT ()
		.test.y.ok (y)
		.iterate.uv (.qfuv.cat.eval, .$levels, .$.PROBS, u=y)
	}
	.catuv (f, "qfuv", "qfuv.cat", x, w, levels)
}

pmfc.cat.cks = function (x, y, sequence,
	kernel=biweight.kernel, bw, smoothness=0.65, w=NA,
	levels, at)
{	fx = pmfuv.cat (x, sequence, w, levels)
	fy = pdfuv.cks (y, FALSE, 30, kernel, bw, smoothness, w)
	nlevels = fx %$% "nlevels"
	ds = categorical.set (pdfuv.cks, y, FALSE, 30, kernel, bw, smoothness, w, group.by=x)
	at = as.numeric (at)
	if (length (at) != 1)
		stop ("at needs to be numeric scalar")

	probs = numeric (nlevels)
	for (i in 1:nlevels)
		probs [i] = fx (i) * ds [[i]](at) / fy (at)
	PROBS = cumsum (probs)
	PROBS [nlevels] = 1

	fx %$% ".probs" = probs
	fx %$% ".PROBS" = PROBS
	fx
}

.val.cat.x.uv = function (nlevels, levels, x)
{	if (is.integer (x) )
	{	if (any (x < 1 | x > nlevels) )
			stop ("x needs to be in [1, nlevels]")
		x
	}
	else if (is.character (x) )
		match (x, levels)
	else if (is.numeric (x) )
		.val.cat.x.uv (nlevels, levels, .int.cat (x ) )
	else
		.cat.err ()
}

.pmfuv.cat.eval.scalar = function (is.weighted, n, x, w, u)
{	if (is.weighted)
		sum (w [u == x])
	else
		sum (u == x) / n
}

.qfuv.cat.eval = function (levels, PROBS, y)
	levels [y <= PROBS][1]
