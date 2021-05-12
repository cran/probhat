#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

ph.mean = function (sf, ..., n.intervals=200) moment (sf, 1, n.intervals=n.intervals)
ph.sd = function (sf, ..., n.intervals=200) sqrt (ph.var (sf, n.intervals=n.intervals) )
ph.var = function (sf, ..., n.intervals=200) moment (sf, 2, n.intervals=n.intervals)
ph.skewness = function (sf, ..., n.intervals=200) moment (sf, 3, n.intervals=n.intervals)
ph.kurtosis = function (sf, ..., n.intervals=200) moment (sf, 4, n.intervals=n.intervals)

moment = function (sf, nth, ..., n.intervals=200)
{	if (nth == 0) 1
	else if (nth == 1) raw.moment (sf, 1, n.intervals=n.intervals)
	else if (nth == 2) central.moment (sf, 2, n.intervals=n.intervals)
	else standardized.moment (sf, nth, n.intervals=n.intervals)
}

central.moment = function (sf, nth, ..., n.intervals=200)
{	if (nth < 1)
		stop ("standardized.moment needs nth >= 1")
	else if (nth == 1)
		0
	else
	{	mean = raw.moment (sf, 1, n.intervals=n.intervals)
		raw.moment (sf, nth, mean, n.intervals=n.intervals)
	}
}

standardized.moment = function (sf, nth, ..., n.intervals=200)
{	if (nth < 2)
		stop ("standardized.moment needs nth >= 2")
	else if (nth == 2)
		1
	else
	{	mean = raw.moment (sf, 1, n.intervals=n.intervals)
		var = raw.moment (sf, 2, mean, n.intervals=n.intervals)
		cen.mmnt = raw.moment (sf, nth, mean, n.intervals=n.intervals)
		cen.mmnt / var^(nth / 2)
	}
}	

raw.moment = function (sf, nth, about=0, ..., n.intervals=200)
{	if (is.pmfuv (sf) )
	{	x = seq (sf)
		y = sf (x)
	}
	else if (is.ccdfuv (sf) )
	{	x = seq (sf, n = n.intervals + 1)
		y = diff (sf (x) )
		x = x [1:n.intervals] + (x [2] - x [1]) / 2
	}
	else
		stop ("needs uv PMF/CCDF")
	if (about != 0)
			x = x - about
	if (nth != 1)
			x = x^nth
	sum (x * y)
}

ntile.names = function (n, symbol="q", ..., emph = n / 2)
{	str = c ("min", paste0 (symbol, 1:(n - 1) ), "max")
	I = unique (c (floor (emph), ceiling (emph) ) )
	str [I + 1] = paste0 ("(", symbol, I, ")")
	str
}

quartiles = function (xf, col=FALSE, ..., prob=FALSE, names = ntile.names (4, "Q", emph=emph), emph=2)
	ntiles (4, xf, col, ..., prob=prob, names=names)
deciles = function (xf, col=FALSE, ..., prob=FALSE, names = ntile.names (10, "D", emph=emph), emph=5)
	ntiles (10, xf, col, ..., prob=prob, names=names)

ntiles = function (n, xf, col=FALSE, ..., prob=FALSE, names)
{	sf = .xqf (xf, ...)

	p = seq (0, 1, length.out = n + 1)
	if (prob)
		names = round (p, 4)
	else
	{	if (missing (names) )
			names = ntile.names (n)
	}
	q = sf (p)
	names (q) = names
	if (col)
	{	q = cbind (q)
		colnames (q) = NULL
	}
	q
}

ph.median = function (xf, ...) ph.quantile (xf, 0.5, ...)

ph.quantile = function (xf, p, ...)
{	sf = .xqf (xf, ...)
	sf (p)
}

iqr = function (xf, P=0.5, ...)
{	sf = .xqf (xf, ...)
	a = (1 - P) / 2
	b = 1 - a

	sf (b) - sf (a)
}

ph.mode = function (sf, infv=FALSE, ..., level.names=FALSE, freq=FALSE, n)
{	if (is.pmfuv (sf) )
	{	.probs = sf %$% ".probs"
		I = which.max (.probs)
		if (infv)
		{	if (freq)
			{	x = seq (sf)[I]
				sf (x, freq=TRUE, n=n)
			}
			else
				.probs [I]
		}
		else
		{	if (is.cat (sf) && level.names)
				(sf %$% "levels")[I]
			else
				seq (sf)[I]
		}
	}
	else if (is.phspline (sf) )
	{	if (freq)
			stop ("freq needs to be false, for continuous models")
		x = ph.modes (sf)
		if (length (x) == 0)
			stop ("no (non-boundary) mode")
		else
		{	y = sf (x)
			I = which.max (y)
			if (infv)
				y [I]
			else
				x [I]
		}
	}
	else
		stop ("needs univariate PMF or spline-based PDF")
}

ph.modes = function (sf, infv=FALSE)
{	if (is.phspline (sf) )
	{	x = argmaxs.chs (sf %$% "spline.function")
		if (infv)
			sf (x)
		else
			x
	}
	else
		stop ("currently, only supports spline-based uv PDFs")
}

rng = function (xf, n=1, ...)
{	if (is.dqf (xf) || is.cqf (xf) )
		y = runif (n)
	else if (is.cchqf (xf) )
		y = matrix (runif (n * xf %$% "m"), nrow=n)
	else
	{	x = as.numeric (xf)
		xf = qfuv.el (x, ...)
		y = runif (n)
	}
	xf (y)
}

.xqf = function (xf, ...)
{	if (is.function (xf) )
	{	if (is.dqf (xf) || is.cqf (xf) )
			xf
		else
			stop ("needs numeric vector or quantile function")
	}
	else
	{	x = as.numeric (xf)
		qfuv.el (x, ...)
	}
}
