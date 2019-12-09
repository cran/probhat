#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.plot.dpd = function (n, xb, h, line.width, line.color, area.color, hmin=0)
{	if (missing (n) )
		n = length (h)
	if (missing (xb) )
		xb = 0:n
	N = 2 * n + 2
	u = v = rep (hmin, N)
	u [1] = xb [1]
	for (i in 1:n)
	{	I = 2 * i
		u [I] = xb [i]
		u [I + 1] = xb [i + 1]
		v [I] = v [I + 1] = h [i]
	}
	u [N] = xb [n + 1]
	polygon (u, v, lwd=line.width, border=line.color, col=area.color)
}

.plot.bars = function (n, xb, h, line.width, line.color, area.color, hmin=0, space=0)
{	if (missing (n) )
		n = length (h)
	if (missing (xb) )
		xb = 0:n
	for (i in 1:n)
	{	u = c (xb [i] + space, xb [i] + space, xb [i + 1] - space, xb [i + 1] - space)
		v = c (hmin, h [i], h [i], hmin)
		polygon (u, v, lwd=line.width, border=line.color, col=area.color)
	}
}

plot.dkernel = function (x,
	main, xlab="x", ylab,
	line.width=1, line.color="black", area.color, 
	xlim, ylim, ..., cdf=FALSE, bw=13)
{	kernel = x
	
	if (missing (main) )
		main = kernel$name
	if (missing (ylab) )
	{	if (cdf)
			ylab = "probability"
		else
			ylab = "mass"
	}
	bw = as.integer (bw)
	if (bw %% 2 == 0)
		bw = bw + 1
	hbw = (bw - 1) / 2
	if (missing (xlim) )
		xlim = c (-hbw, hbw)
	x = xlim [1]:xlim [2]
	if (cdf)
		y = kernel$cdf (bw, x)
	else
		y = kernel$pmf (bw, x)
	xlim [1] = xlim [1] - 1
	if (missing (ylim) )
		ylim = c (0, 1.025 * max (y) )
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color
	}
	plot.new ()
	plot.window (xlim=xlim, ylim=ylim, yaxs="i")
	title (main,, xlab, ylab)
	.plot.dpd (,c (x [1] - 1, x), y, line.width, line.color, area.color)
	box ()
	axis (1)
	axis (2)
}

.plot.dksuv = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color, ..., axis.ticks=TRUE, xlim, ylim)
{	x = seq (f)
	y = f (x)
	if (missing (main) )
		main = ""
	if (missing (xlim) )
	{	xlim = range (f)
		xlim [1] = xlim [1] - 1
	}
	if (missing (ylim) )
		ylim = c (0, 1.025 * max (y) )
	plot.new ()
	plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
	title (main,, xlab, ylab)
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color
	}
	.plot.dpd (,c (x [1] - 1, x), y, line.width, line.color, area.color)
	box ()
	if (axis.ticks)
	{	axis (1, x - 0.5, .pblabs (x, f %$% x) )
		axis (2)
	}
}

.plot.dksuv.qf = function (F.inv,
	main, xlab, ylab,
	line.width, line.color, area.color, ...)
{	if (missing (main) )
		main = ""
	xlim = F.inv %$% "xlim"
	n = diff (xlim) + 1
	x = xlim [1]:xlim [2]
	y = F.inv %$% ".PROBS"
	yb = c (0, y)
	xlim [1] = xlim [1] - 1
	xlim [2] = xlim [2] + 0.025 * diff (xlim)
	plot.new ()
	plot.window (xlim = c (0, 1), ylim=xlim, xaxs="i", yaxs="i")
	title (main,, xlab, ylab)
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color.2
	}
	.plot.dpd (n, yb, x, line.width, line.color, area.color, xlim [1])
	box ()
	axis (1)
	axis (2, x - 0.5, .pblabs (x, F.inv %$% x) )
}

.plot.dksuv.2 = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color)
{	x = seq (f)
	y = f (x)
	if (missing (main) )
		main = ""
	xlim = range (f)
	xlim [1] = xlim [1] - 1
	height.main = max (y)
	height.sub = height.main / 4.5
	space = 0.025 * height.main
	ymin = -1 * height.sub - space
	ymax = height.main + space
	ylim = c (ymin, ymax)
	plot.new ()
	plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
	title (main,, xlab, ylab)
	options = getOption ("probhat")
	if (missing (area.color) )
		area.color = options$area.color
	area.color.2 = options$area.color.2
	.plot.dpd (,c (x [1] - 1, x), y, line.width, line.color, area.color)
	p = .pobs (x, f %$% "x", f %$% "h", height.sub)  + ymin
	.plot.bars (,c (x [1] - 1, x), p, line.width, line.color, area.color.2, ymin)
	abline (h=0)
	y = pretty (y)
	box ()
	axis (1, x - 0.5, .pblabs (x, f %$% "x") )
	axis (2, y [y >= 0])
}

.pobs = function (u, x, h, height.sub)
{	y = rep (0, length (u) )
	y [match (x, u)] = h
	height.sub / max (h) * y
}

.pblabs = function (u, x)
{	rnames = names (x)
	if (is.null (rnames) )
		u
	else
	{	blabs = rep ("", length (u) )
		blabs [match (x, u)] = rnames
		blabs
	}
}

.plot.categorical = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color)
{	if (missing (main) )
		main = ""
	if (missing (xlab) )
		xlab = f %$% "variable.name"
	if (missing (ylab) )
		ylab = "probability"
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color
	}
	at = 1:(f %$% "nlevels")
	x = f %$% "levels"
	y = f (x)
	ylim = c (0, 1.025 * max (y) )
	plot.new ()
	plot.window (xlim = c (0, f %$% "nlevels"), ylim=ylim, yaxs="i")
	.plot.bars (,,y, line.width, line.color, area.color,,space=0.05)
	title (main,, xlab, ylab)
	box ()
	axis (1, at - 0.5, x)
	axis (2)
}

.plot.categorical.qf = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color)
{	if (missing (main) )
		main = ""
	if (missing (xlab) )
		xlab = "probability"
	if (missing (ylab) )
		ylab = f %$% "variable.name"
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color.2
	}
	n = f %$% "nlevels"
	at = 1:n
	x = f %$% "levels"
	yb = c (0, f %$% ".PROBS")
	plot.new ()
	plot.window (xlim = c (0, 1), ylim = c (0, 1.025 * n), yaxs="i")
	.plot.bars (n, yb, at, line.width, line.color, area.color,, space=0.05)
	title (main,, xlab, ylab)
	box ()
	axis (1)
	axis (2, at - 0.5, x)
}
