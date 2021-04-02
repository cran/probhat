#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.dksuv = function (f, classes, bw.method, kernel, bw, smoothness, XLIM, x, h, tail)
{	xname = .varname (x)
	blabs = .blabs (x)
	x = .val.x.uv (x, TRUE)
	n = length (x)
	names (x) = blabs
	xlim = range (x)
	dx = diff (x)

	h = .val.hvec (n, h)

	I = order (x)
	x = x [I]
	h = h [I]

	ux = unique (x)
	nu = length (ux)
	if (nu != n)
	{	n = nu
		h = .aggregate (nu, ux, x, h)
		x = ux
	} 
	.xsum = sum (h)

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
	hbw = (bw - 1) / 2
	kernel = .val.k (kernel)
	kernel = discretized.kernel (bw, kernel)

	XLIM = .val.XLIM.uv (XLIM, TRUE)
	.is.trunc = .is.trunc.uv (hbw, XLIM, xlim)
	.any.trunc = any (.is.trunc)

	if (.any.trunc)
	{	.check.x.inside.uv (.is.trunc, XLIM, x)

		xlim = .update.xlim (.is.trunc, XLIM, xlim, TRUE)
		fdata = .update.x.uv (.is.trunc, hbw, XLIM, n, x, h, TRUE, TRUE)
	}
	else
		fdata = list (n=n, x=x, h=h)
	w = fdata$h / sum (fdata$h)

	u = xlim [1]:xlim [2]
	.probs = .iterate.uv (.pmfuv.dks.eval.0, kernel@f, fdata$x, w, u=u)
	.probs = .probs / sum (.probs)

	if (missing (tail) )
	{	.low = tail = NA
		.PROBS = .cumsum2 (.probs, FALSE)
	}
	else
	{	tail = .val.tail (tail)
		.low = (tail == "lower")
		.PROBS = .cumsum2 (.probs, ! .low)
	}
	
	.EXTEND (f, classes,
		.any.trunc, .is.trunc,
		.probs, .PROBS, .low, .xsum,
		xname, XLIM, kernel, tail, bw, smoothness,
		n, xlim, x, h)
}

pmfuv.dks = function (x = 1:length (h), h=1, ...,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default",
	Xlim = c (a, b),
	a = min1 (x), b=Inf)
{	.arg.error (...)
	.dksuv (.pmfuv.dks.eval, .CV.pmfuv.dks, bw.method, kernel, bw, smoothness, Xlim, x, h)
}

cdfuv.dks = function (x = 1:length (h), h=1, ...,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", tail="lower",
	Xlim = c (a, b),
	a = min1 (x), b=Inf)
{	.arg.error (...)
	.dksuv (.cdfuv.dks.eval, .CV.cdfuv.dks, bw.method, kernel, bw, smoothness, Xlim, x, h, tail)
}

qfuv.dks = function (x = 1:length (h), h=1, ...,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default",
	Xlim = c (a, b),
	a = min1 (x), b=Inf)
{	.arg.error (...)
	.dksuv (.qfuv.dks.eval, .CV.qfuv.dks, bw.method, kernel, bw, smoothness, Xlim, x, h, tail="lower")
}

min1 = function (x)
{	minx = min (x)
	if (minx == 1) 1
	else 0
}

.aggregate = function (nu, ux, ox, oh)
{	h = numeric (nu)
	for (i in 1:nu)
	{	I = (ux [i] == ox)
		h [i] = sum (oh [I])
	}
	h
}
