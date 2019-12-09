#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.plot.cpd = function (n, x, y, line.width, line.color, area.color, ymin=0)
{	x = c (x [1], x, x [n])
	y = c (ymin, y, ymin)
	polygon (x, y, lwd=line.width, border=line.color, col=area.color)
}

.empty.kernel.plot = function (main="", y)
{	plot.new ()
	plot.window (xlim = c (-1.05, 1.05), ylim = c (0, 1.25), xaxs="i", yaxs="i")
	title (main)
	abline (h=y, lty=3)
}

plot.ckernel = function (x,
	main, xlab="x", ylab,
	line.width=1, line.color="black", area.color, 
	xlim = c (-1.05, 1.05), ylim,
	..., add=FALSE, cdf=FALSE)
{	kernel = x

	options = getOption ("probhat")
	if (missing (main) )
		main = kernel$name
	x = seq (-1, 1, length.out=200)
	if (cdf)
	{	if (missing (ylab) )
			ylab = "probability"
		y = kernel$cdf (x)
	}
	else
	{	if (missing (ylab) )
			ylab = "density"
		y = kernel$pdf (x)
	}
	if (missing (area.color) )
		area.color = options$area.color
	if (! add)
	{	if (missing (ylim) )
			ylim = c (0, 1.025 * max (y) )
		plot.new ()
		plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
		title (main,, xlab, ylab)
	}
	.plot.cpd (200, x, y, line.width, line.color, area.color)
	box ()
	if (! add)
	{	axis (1)
		axis (2)
	}
}

.plot.cksuv = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color, ..., is.qf=FALSE, axis.ticks=TRUE, xlim, ylim)
{	if (is.qf)
		x = seq (0, 1, length.out=200)
	else
		x = seq (f, 200)
	y = f (x)
	if (missing (main) )
		main = ""
	if (missing (xlim) )
	{	if (is.qf)
			xlim = c (0, 1)
		else
			xlim = range (f)
	}
	if (missing (ylim) )
	{	if (is.qf)
		{	ylim = range (y)
			os = 0.025 * diff (ylim)
			ylim = c (ylim [1], ylim [2] + os)
		}
		else
			ylim = c (0, 1.025 * max (y) )
	}
	plot.new ()
	plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
	title (main,, xlab, ylab)
	if (missing (area.color) )
	{	options = getOption ("probhat")
		if (is.qf)
			area.color = options$area.color.2
		else
			area.color = options$area.color
	}
	if (is.qf)
		ymin = min (y)
	else
		ymin = 0
	.plot.cpd (200, x, y, line.width, line.color, area.color, ymin)
	box ()
	if (axis.ticks)
	{	axis (1)
		axis (2)
	}
}

.add.plot.cksuv = function (f, line.width, line.color, area.color)
{	x = seq (f, 200)
	y = f (x)
	.plot.cpd (200, x, y, line.width, line.color, area.color)
}

.plot.cksuv.2 = function (f,
	main, xlab, ylab,
	line.width, line.color, area.color)
{	xlim = range (f)
	x = seq (f, 200)
	y = f (x)
	if (missing (main) )
		main = ""
	top = ymax = max (y)
	space = 0.025 * max (y)
	ymin = ymax / (-4.5)
	ymax = ymax + space
	if (missing (area.color) )
	{	options = getOption ("probhat")
		area.color = options$area.color
	}
	plot.new ()
	plot.window (xlim=xlim, ylim = c (ymin, ymax), xaxs="i", yaxs="i")
	title (main,, xlab, ylab)
	.plot.cpd (200, x, y, line.width, line.color, area.color)

	is.weighted = f %$% is.weighted
	n = f %$% n
	x = f %$% x
	w = f %$% w
	y = runif (n, ymin + space, - space)
	.plot.points.2 (is.weighted, x, y, w)

	abline (h=0)
	box ()
	axis (1)
	axis (2, round (seq (0, top, length.out=4), 2) )
}

.plot.points = function (is.weighted, x, y, w, xlab, ylab)
{	plot (x, y, xlab=xlab, ylab=ylab, pch=NA)
	.plot.points.2 (is.weighted, x, y, w)
}

.plot.points.2 = function (is.weighted, x, y, w)
{	if (is.weighted)
	{	v = 1 - w
		points (x, y, col = rgb (v, v, v) )
	}
	else
	{	options = getOption ("probhat")
		col=options$area.color.2
		points (x, y, pch=16, col=col)
	}
	points (x, y)
}
