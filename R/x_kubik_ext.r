#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.intset = function (x)
{	dx = diff (x)
	dx [dx < 1e-6] = 0
	I = (dx > 0)
	k1 = diff (c (0, I) )
	k2 = diff (c (I, 0) )
	cbind (which (k1 == 1), which (k2 == -1) + 1)
}

.incr.chs = function (cx, cy, cb, ...)
{	constraints = chs.constraints (increasing=TRUE)
	chs (cx, cy, cb, constraints=constraints, ...)
}

.spline = function (f, with.cdf, nc)
{	cx = seq (f, n=nc)
	cy = f (cx)
	if (with.cdf)
	{	cy [c (1, nc)] = c (0, 1)
		.incr.chs (cx, cy, outside = c (0, 1) )
	}
	else
	{	cy [c (1, nc)] = c (0, 0)
		chs (cx, cy, outside = c (0, 0) )
	}	
}

.modified.spline.transposed = function (cx, cy)
{	nc = length (cx)
	cy [1] = 0
	cy [nc] = 1

	ints = .intset (cy)
	nsplines = nrow (ints)
	if (nsplines == 1)
	{	I = ints [1, 1]:ints [1, 2]
		.incr.chs (cy [I], cx [I])
	}
	else
	{	knots = cy [ints [,1][-1]]
		splines = vector ("list", nsplines)
		for (i in 1:nsplines)
		{	I = ints [i, 1]:ints [i, 2]
			splines [[i]] = .incr.chs (cy [I], cx [I])
		}
		f = function (x)
		{	. = THAT ()
			.iterate.uv (.nested.chs.eval, .$knots, .$splines, u=x)
		}
		EXTEND (f, "nested.chs", nsplines, knots, splines)
	}
}

.nested.chs.eval = function (knots, splines, x)
{	I = 1 + sum (x >= knots)
	splines [[I]](x)
}
