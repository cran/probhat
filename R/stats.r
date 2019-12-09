#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

ph.mean = function (f, n=200)
	ph.moment (f, 1, n)

ph.sd = function (f, n=200)
	sqrt (ph.var (f, n) )

ph.var = function (f, n=200)
	ph.moment (f, 2, n)

ph.skewness = function (f, n=200)
	ph.moment (f, 3, n)

ph.kurtosis = function (f, n=200)
	ph.moment (f, 4, n)

ph.median = function (F.inv)
	ph.quantile (F.inv, 0.5)

ph.quantile = function (F.inv, p)
{	if (is.qfuv (F.inv) )
		F.inv (p)
	else
		stop ("needs univariate quantile function")
}

ph.mode = function (f)
{	if (is.pmfuv (f) )
		(f %$% "x")[which.max (f %$% ".probs")]
	else if ( (inherits (f, "pdfuv.cks") || inherits (f, "pdfc.cks") ) && f %$% "is.spline")
	{	x = chs.argmaxs (f %$% "spline.function")
		if (length (x) == 0)
			x = f %$% "xlim"
		y = f (x)
		x [which.max (y)]
	}
	else
		stop ("needs univariate pmf or pdf")
}

ph.moment = function (f, nth, n=200)
{	if (nth == 0)
		1
	else
	{	mean = ph.moment.2 (f, 1, 0, n)
		if (nth == 1)
			mean
		else
		{	var = ph.moment.2 (f, 2, mean, n)
			if (nth == 2)
				var
			else
			{	um = ph.moment.2 (f, nth, mean, n)
				um / var ^ (nth / 2)
			}
		}
	}
}

ph.moment.2 = function (f, nth=1, about=0, n=200)
{	if (is.pmfuv (f) )
	{	x = seq (f)
		y = f %$% ".probs"
	}
	else if (is.cpd (f) && is.cdfuv (f) )
	{	x = seq (f, n + 1)
		y = diff (f (x) )
		x = x [1:n] + (x [2] - x [1]) / 2
	}
	else
		stop ("needs univariate pmf or continuous cdf")
	if (about != 0)
			x = x - about
	if (nth != 1)
			x = x^nth
	sum (x * y)
}

ph.rng = function (F.inv, n=1)
{	if (is.qfuv (F.inv) )
		y = runif (n)
	else if (is.chqf (F.inv) )
		y = matrix (runif (n * F.inv %$% "m"), nrow=n)
	else
		stop ("needs univariate qf, or chqf.cks object")
	F.inv (y)
}
