#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.val.XLIM.uv = function (XLIM, as.int=FALSE)
{	XLIM = as.vector (XLIM)
	if (length (XLIM) != 2)
		stop ("vector limits should be length two")
	if (! is.numeric (XLIM) )
		stop ("limits need to be integer/numeric")
	if (XLIM [2] > XLIM [1])
	{	if (as.int)
		{	I = (is.finite (XLIM) )
			if (sum (I) > 0)
			{	sub = XLIM [I]
				rsub = round (sub)
				if (rsub != sub)
					stop ("discrete limits, needs integer-equivalent values")
			}
		}
		XLIM
	}
	else
		stop ("lower limits >= upper limits")
	
}

.val.XLIM.mv = function (XLIM, m)
{	if (! is.numeric (XLIM) )
		stop ("limits need to be numeric")
	if (is.matrix (XLIM) )
	{	if (ncol (XLIM) != 2)
			stop ("matrix limits need two columns")
		nr = nrow (XLIM)
		if (nr == 1)
			XLIM = cbind (rep (XLIM [1, 1], m), rep (XLIM [1, 2], m) )
		else if (nr != m)
			stop ("matrix limits need 1 or m rows")
		if (all (XLIM [,2] > XLIM [,1]) )
			XLIM
		else
			stop ("lower limits >= upper limits")
	}
	else
	{	XLIM = as.vector (XLIM)
		if (length (XLIM) != 2)
			stop ("vector limits should be length two")
		.val.XLIM.mv (rbind (XLIM), m)
	}
}

.is.trunc.uv = function (hbw, XLIM, xlim0)
{	lo = (xlim0 [1] - hbw < XLIM [1])
	up = (xlim0 [2] + hbw > XLIM [2])
	c (lo, up)
}

.is.trunc.mv = function (hbw, XLIM, xlim0)
{	lo = (xlim0 [,1] - hbw < XLIM [,1])
	up = (xlim0 [,2] + hbw > XLIM [,2])
	cbind (lo, up)
}

.update.xlim = function (ist, XLIM, xlim, as.int=FALSE)
{	xlim [ist] = XLIM [ist]
	if (as.int)
		xlim = as.integer (xlim)
	xlim
}

.update.x.uv = function (ist, hbw, XLIM, n, x, w=NA, isw, as.int=FALSE)
{	n.lo = n.up = 0
	x.lo = x.up = w.lo = w.up = numeric ()
	if (ist [1])
	{	I.lo = (
			x != XLIM [1] &
			x < XLIM [1] + hbw)
		n.lo = sum (I.lo)
		x.lo = 2 * XLIM [1] - x [I.lo]
		if (isw)
			w.lo = w [I.lo]
	}
	if (ist [2])
	{	I.up = (
			x != XLIM [2] &
			x > XLIM [2] - hbw)
		n.up = sum (I.up)
		x.up = 2 * XLIM [2] - x [I.up]
		if (isw)
			w.up = w [I.up]
	}
	n = n + n.lo + n.up
	x = c (x, x.lo, x.up)
	if (isw)
	{	w = c (w, w.lo, w.up)
		w = w / sum (w)
	}
	if (as.int)
		list (n=n, x = as.integer (x), h=w)
	else
		list (n=n, x=x, w=w)
}

.update.x.mv = function (ist, hbw, XLIM, n, m, x, w, isw, as.int=FALSE)
{	IST = apply (ist, 2, any)

	n.lo = n.up = 0
	x.lo = x.up = w.lo = w.up = numeric ()
	if (IST [1])
	{	I.lo = rep (FALSE, n)
		for (j in which (ist [,1]) )
		{	I.lo = I.lo | (
				x [,j] != XLIM [j, 1] &
				x [,j] < XLIM [j, 1] + hbw [j])
		}
		n.lo = sum (I.lo)
		x.lo = .reflect.mv (ist [,1], XLIM [,1], 1, n.lo, x [I.lo,, drop=FALSE])
		if (isw)
			w.lo = w [I.lo]
	}
	if (IST [2])
	{	I.up = rep (FALSE, n)
		for (j in which (ist [,2]) )
		{	I.up = I.up | (
				x [,j] != XLIM [j, 2] &
				x [,j] > XLIM [j, 2] - hbw [j])
		}
		n.up = sum (I.up)
		x.up = .reflect.mv (ist [,2], XLIM [,2], -1, n.up, x [I.up,, drop=FALSE])
		if (isw)
			w.up = w [I.up]
	}

	n = n + n.lo + n.up
	x = rbind (x, x.lo, x.up)
	if (isw)
	{	w = c (w, w.lo, w.up)
		w = w / sum (w)
	}
	list (n=n, x=x, w=w)
}

#all points reflected, but not necessarily all variables
.reflect.mv = function (ist1, XLIM1, dvec, n, x)
{	nref = sum (ist1)

	offset = matrix (rep (XLIM1 [ist1], each=n), n, nref)
	xsub = x [,ist1, drop=FALSE]

	dist = xsub - offset
	mins = apply (dvec * dist, 1, min)
	mins = matrix (rep (mins, times=nref), n, nref)
	p0 = offset + dist - dvec * mins

	x [,ist1] = 2 * p0 - xsub
	x
}

.sub.data.col = function (data, j)
{	data$xlim = data$xlim [j,, drop=FALSE]
	data$x = data$x [,j, drop=FALSE]
	data
}

.check.x.inside.uv = function (ist, XLIM, x, str="x")
	.check.x.inside.mv (rbind (ist), rbind (XLIM), 1, cbind (x), str)

.check.x.inside.mv = function (ist, XLIM, m, x, str="x")
{	for (j in 1:m)
	{	if (ist [j, 1] && any (x [,j] < XLIM [j, 1]) )
			stop (sprintf ("%s values below lower limit", str) )
		if (ist [j, 2] && any (x [,j] > XLIM [j, 2]) )
			stop (sprintf ("%s values above upper limit", str) )
	}
}

.val.u.uv = function (u, anytr=FALSE, ist, XLIM)
{	if (anytr)
		.check.x.inside.uv (ist, XLIM, u, "eval point")
	u
}

.val.u.mv = function (m, u, anytr=FALSE, ist, XLIM)
{	if (! is.matrix (u) )
		u = rbind (u)
	if (m != ncol (u) )
		stop ("incorrect number of columns")
	if (anytr)
		.check.x.inside.mv (ist, XLIM, m, u, "eval point")
	u
}

.is.trunc.lower.side = function (ist, low)
	( (ist [,1] & low) | (ist [,2] & ! low) )

.cdf.lower.side = function (low, xlim)
{	u = xlim [,2]
	u [low] = xlim [low ,1]
	u
}

.cdf.upper.side = function (low, xlim)
{	u = xlim [,2]
	u [! low] = xlim [! low,1]
	u
}

.cdfv.lower.side = function (Fh, low, isltr, M, xlim, ..., hi)
{	lo = .cdf.lower.side (low, xlim)
	if (missing (hi) )
		hi = .cdf.upper.side (low, xlim)
	ntr = sum (isltr)
	u = matrix (rep (lo, each=ntr), ntr, M)
	J = which (isltr)
	for (j in 1:ntr)
	{	k = J [j]
		u [j, k] = hi [k]
	}
	max (Fh (u, .ignore.trunc=TRUE) )
}

.cdf.scalef = function (Fh, M=1,
	low = Fh %$% ".low",
	anyltr = Fh %$% ".any.trunc.lower",
	isltr = Fh %$% ".is.trunc.lower",
	xlim = (Fh %$% "data")$xlim)
{	isltr = rbind (isltr)
	xlim = rbind (xlim)
	nr = nrow (xlim)
	J = (1 + nr - M):nr
	xlim = xlim [J,, drop=FALSE]
	hi = .cdf.upper.side (low, xlim)
	fhi = Fh (rbind (hi), .ignore.trunc=TRUE) 
	if (anyltr)
		flo = .cdfv.lower.side (Fh, low, isltr, M, xlim, hi=hi)
	else
		flo = 0
	1 / (fhi - flo)
}
