#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.el = function (f, classes, x, w, inverse)
{	xname = .varname (x)
	x = .val.x.uv (x)
	n = length (x)
	if (n < 2)
		stop ("need 2 or more obs")
	randomize = (.any.duplicates (x) )
	if (randomize)
	{	sd = sd (x) / 1000
		while (randomize)
		{	x = x + runif (n, -sd, sd)
			randomize = (.any.duplicates (x) )
		}
	}
	is.weighted = (! (missing (w) || is.na (w [1]) ) )
	x.order = order (x)
	x = x [x.order]
	xlim = c (x [1], x [n])

	if (is.weighted)
	{	w = .val.w (is.weighted, n, w)
		w = w [x.order]
		intw = (w [-n] + w [-1]) / 2
		outw = (w [1] + w [n]) / 2
		intw =  intw / (1 - outw)
		y = cumsum (c (0, intw) )
	}
	else
	{	w = NA
		y = (0:(n - 1) ) / (n - 1)
	}
	y [n] = 1
	if (inverse)
		spline.function = .incr.chs (y, x)
	else
		spline.function = .incr.chs (x, y, outside = c (0, 1) )
	.EXTEND (f, classes,
		xname, is.weighted, spline.function, xlim, n, x, w)
}
 
cdfuv.el = function (x, ..., w)
{	F = function (x)
	{	. = .THAT ()
		.$spline.function (x)
	}
	.arg.error (...)
	.el (F, .CV.cdfuv.el, x, w, FALSE)
}

qfuv.el = function (x, ..., w)
{	F.inv = function (p)
	{	. = .THAT ()
		.test.y.ok (p)
		.$spline.function (p)
	}
	.arg.error (...)
  	.el (F.inv, .CV.qfuv.el, x, w, TRUE)
}
