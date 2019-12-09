#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.el = function (f, superclass, subclass, x, w, inverse=FALSE)
{	variable.name = .varname (x)
	x = .val.x.uv (x)
	n = length (x)
	randomize = (.any.duplicates (x) )
	if (randomize)
	{	sd = sd (x) / 1000
		while (randomize)
		{	x = x + runif (n, -sd, sd)
			randomize = (.any.duplicates (x) )
		}
	}
	is.weighted = (!is.na (w [1]) )
	x.order = order (x)
	x = x [x.order]
	xlim = c (x [1], x [n])
	if (is.weighted)
	{	w = .val.w (n, w)
		w = w [x.order]
		intw = (w [-n] + w [-1]) / 2
		outw = (w [1] + w [n]) / 2
		intw =  intw / (1 - outw)
		y = cumsum (c (0, intw) )
	}
	else
		y = (0:(n - 1) ) / (n - 1)
	y [n] = 1
	if (inverse)
		spline.function = chs (y, x)
	else
		spline.function = chs (x, y, outside = c (0, 1) )
	EXTEND (f, c (subclass, "el", superclass, "cpd", "phmodel"),
		variable.name, is.weighted, spline.function, xlim, n, x, w)
}
 
cdf.el = function (x, w=NA)
{	F = function (x)
	{	. = THAT ()
		.$spline.function (x)
	}
	.el (F, "cdfuv", "cdf.el", x, w)
}

qf.el = function (x, w=NA)
{	F.inv = function (y)
	{	. = THAT ()
		.test.y.ok (y)
		.$spline.function (y)
	}
  	.el (F.inv, "qfuv", "qf.el", x, w, TRUE)
}
