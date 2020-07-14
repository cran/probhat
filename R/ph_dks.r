#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.dksuv = function (f, classes, bw.method, kernel, bounded, lower, upper, bw, smoothness, x, h, freq=FALSE)
{	if (missing (x) )
	{	if (missing (h) )
			stop ("either x or h needed")
		else
		{	variable.name = "x"
			h = as.numeric (h)
			n = length (h)
			if (n == 0)
				stop ("x missing and length (h) == 0")
			xlim = c (1, n)
			x = 1:n

			warning ("x missing, default 1:length (h)")
		}
	}
	else
	{	variable.name = .varname (x)
		blabs = .blabs (x)
		x = .val.x.uv (x, TRUE)
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
	N = sum (h)
	if (any (h < 0) )
		stop ("h needs non-negative values")

	xrng = diff (xlim) + 1
	if (missing (bw) )
	{	bw = auto.dbw (rep (x, times=h), bw.method=bw.method, smoothness=smoothness)
		bw = as.integer (round (bw) )
		if (bw %% 2 == 0)
			bw = bw + 1
	}
	else
	{	if (bw < 1)
			stop ("bw needs to be positive integer")
		if (bw %% 2 == 0)
		{	warning ("bw even, so incremented")
			bw = bw + 1
		}
		smoothness = NA
	}
	kernel = .val.k (kernel, bw, super="dkernel")
	hbw = (bw - 1) / 2

	bounded = rep_len (bounded, 2)
	if (bounded [1])
	{	if (! missing (lower) )
		{	lower = as.integer (lower)
			if (xlim [1] < lower) stop ("min (x) < lower")
			else	xlim [1] = lower
		}
	}
	else
		xlim [1] = xlim [1] - hbw
	if (bounded [2])
	{	if (! missing (upper) )
		{	upper = as.integer (upper)
			if (xlim [2] > upper) stop ("max (x) > upper")
			else	xlim [2] = upper
		}
	}
	else
		xlim [2] = xlim [2] + hbw

	u = xlim [1]:xlim [2]
	if (any (bounded) )
	{	rx = rev (x)
		rh = rev (h)
		x2 = c (2 * xlim [1] - rx [-n], x, 2 * xlim [2] - rx [-1])
		h2 = c (rh [-n], h, rh [-1])
		.probs = .iterate.uv (.pmfuv.dks.eval.2, kernel$pmf, x2, h2, u=u)
	}
	else
		.probs = .iterate.uv (.pmfuv.dks.eval.2, kernel$pmf, x, h, u=u)
	.probs = .probs / sum (.probs)

	.PROBS = cumsum (.probs)
	.PROBS [length (u)] = 1

	EXTEND (f, classes,
		.probs, .PROBS, N,
		variable.name, kernel, bounded, freq, bw, smoothness,
		xlim, n, x, h)
}

.pmfuv.dks.f = function (., x, freq)
{	x = .val.x.uv (x)
	y = .iterate.uv (.pmfuv.dks.eval, .$xlim, .$.probs, u=x)
	if (freq)
		.$N * y
	else
		y
}

pmfuv.dks = function (x, h, ..., bw.method="ph.default",
	kernel=binomial.dkernel, bounded = c (TRUE, FALSE), freq=FALSE, lower, upper, bw, smoothness=1)
{	f = function (x, ..., freq)
	{	. = THAT ()
		if (missing (freq) )
			freq = .$freq
		.pmfuv.dks.f (., x, freq)
	}
	.dksuv (f, .CV.pmfuv.dks, bw.method, kernel, bounded, lower, upper, bw, smoothness, x, h, freq)
}

cdfuv.dks = function (x, h, ..., bw.method="ph.default",
	kernel=binomial.dkernel, bounded = c (TRUE, FALSE), lower, upper, bw, smoothness=1)
{	f = function (q)
	{	. = THAT ()
		q = .val.x.uv (q)
		.iterate.uv (.cdfuv.dks.eval, .$xlim, .$.PROBS, u=q)
	}
	.dksuv (f, .CV.cdfuv.dks, bw.method, kernel, bounded, lower, upper, bw, smoothness, x, h)
}

qfuv.dks = function (x, h, ..., bw.method="ph.default",
	kernel=binomial.dkernel, bounded = c (TRUE, FALSE), lower, upper, bw, smoothness=1)
{	f = function (p)
	{	. = THAT ()
		.test.y.ok (p)
		.iterate.uv (.qfuv.dks.eval, .$xlim, .$.PROBS, u=p)
	}
	.dksuv (f, .CV.qfuv.dks, bw.method, kernel, bounded, lower, upper, bw, smoothness, x, h)
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

.pmfuv.dks.eval.2 = function (kpmf, x, w, u)
	sum (w * kpmf (u - x) )

.aggregate = function (nu, ux, ox, oh)
{	h = numeric (nu)
	for (i in 1:nu)
	{	I = (ux [i] == ox)
		h [i] = sum (oh [I])
	}
	h
}
