#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.dksuv = function (f, superclass, subclass, kernel, bw, lower, x, h)
{	if (missing (x) )
	{	if (missing (h) )
			stop ("either x or h needed")
		else
		{	variable.name = "x"
			h = as.numeric (h)
			n = length (h)
			xlim = c (1, n)
			x = 1:n
		}
	}
	else
	{	variable.name = .varname (x)
		blabs = .blabs (x)
		x = .val.x.uv (x)
		n = length (x)
		names (x) = blabs
		xlim = range (x)
		dx = diff (x)
		if (missing (h) )
			h = rep (1, n)
		else
			h = as.integer (h)
		ux = unique (x)
		nu = length (ux)
		if (nu != n)
		{	n = nu
			h = .aggregate (nu, ux, x, h)
			x = ux
		} 
	}
	if (any (h < 0) )
		stop ("h needs non-negative values")
	bw = as.integer (bw)
	if (bw < 1)
		stop ("bw needs to be positive integer")
	if (bw %% 2 == 0)
	{	warning ("bw even, so incremented")
		bw = bw + 1
	}
	kernel = .val.k (kernel)
	hbw = (bw - 1) / 2
	truncated = TRUE
	if (! missing (lower) )
	{	if (lower == -Inf)
		{	truncated = FALSE
			xlim [1] = xlim [1] - hbw
		}
		else
			xlim [1] = lower
	}
	xlim [2] = xlim [2] + hbw
	u = xlim [1]:xlim [2]

	if (truncated)
	{	xlim2 = xlim
		x2 = c (2 * xlim [1] - rev (x [-1]), x)
		h2 = c (rev (h [-1]), h)
		w = h2 / sum (h2)
		truncated.area = .cdfuv.dks.eval.2 (kernel$cdf, bw, xlim2, x2, w)
		.probs = .iterate.uv (.pmfuv.dks.eval.2, kernel$pmf, bw, x2, w, u=u)
		.probs = .probs / (1 - truncated.area)
	}
	else
		.probs = .iterate.uv (.pmfuv.dks.eval.2, kernel$pmf, bw, x, h / sum (h), u=u)
	.PROBS = cumsum (.probs)
	.PROBS [length (u) + 1] = 1

	EXTEND (f, c (subclass, "dksuv", superclass, "dpd", "phmodel"),
		.probs, .PROBS,
		variable.name, kernel, bw,
		xlim, n, x, h)
}

pmfuv.dks = function (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
{	f = function (x)
	{	. = THAT ()
		x = .val.x.uv (x)
		.iterate.uv (.pmfuv.dks.eval, .$xlim, .$.probs, u=x)
	}
	.dksuv (f, "pmfuv", "pmfuv.dks", kernel, bw, lower, x, h)
}

cdfuv.dks = function (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
{	f = function (x)
	{	. = THAT ()
		x = .val.x.uv (x)
		.iterate.uv (.cdfuv.dks.eval, .$xlim, .$.PROBS, u=x)
	}
	.dksuv (f, "cdfuv", "cdfuv.dks", kernel, bw, lower, x, h)
}

qfuv.dks = function (x, h, sequence, kernel=binomial.kernel, bw=1, lower)
{	f = function (y)
	{	. = THAT ()
		.test.y.ok (y)
		.iterate.uv (.qfuv.dks.eval, .$xlim, .$.PROBS, u=y)
	}
	.dksuv (f, "qfuv", "qfuv.dks", kernel, bw, lower, x, h)
}

.pmfuv.dks.eval = function (xlim, probs, x)
{	if (x < xlim [1] || x > xlim [2])
		0
	else
		probs [x - xlim [1] + 1]
}

.cdfuv.dks.eval = function (xlim, PROBS, x)
{	if (x < xlim [1])
		0
	else if (x > xlim [2])
		1
	else
		PROBS [x - xlim [1] + 1]
}

.qfuv.dks.eval = function (xlim, PROBS, y)
	xlim [1] + which (y <= PROBS)[1] - 1

.pmfuv.dks.eval.2 = function (kpmf, bw, x, w, u)
	sum (w * kpmf (bw, u - x) )

.cdfuv.dks.eval.2 = function (kcdf, bw, xlim, x, w)
	sum (w * kcdf (bw, xlim [1] - 1 - x) )

.aggregate = function (nu, ux, ox, oh)
{	h = numeric (nu)
	for (i in 1:nu)
	{	I = (ux [i] == ox)
		h [i] = sum (oh [I])
	}
	h
}
