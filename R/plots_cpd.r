#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.pcpd0 = function (x, y, line.width, line.color, fill.color, ymin=0)
{	x = c (x [1], x, x [length (x)])
	y = c (ymin, y, ymin)
	polygon (x, y, border=NA, col=fill.color)
	lines (x, y, lwd=line.width, col=line.color)
}

plot_cpd = function (sf, data=FALSE, ...,
	main, xlab, ylab,
	xlim, ylim,
	add=FALSE, axes=TRUE,
	line.width, line.color, fill.color)
{	axes = rep_len (axes, 2)
	options = getOption ("probhat")
	if (missing (main) )
		main = ""
	if (missing (line.width) )
		line.width = options$main.line.width
	if (missing (line.color) )
		line.color = options$main.line.color

	if (is.pdf (sf) )
	{	h0 = 0
		if (missing (xlim) )
			xlim = range (sf)
		x = seq (sf)
		y = sf (x)
		if (missing (ylim) )
			ylim = c (0, max (y) )

		if (missing (xlab) )
			xlab = .deflab (sf)
		if (missing (ylab) )
			ylab = "density"

		if (missing (fill.color) )
			fill.color = options$main.fill.color
	}
	else if (is.ccdf (sf) )
	{	h0 = 0
		if (missing (xlim) )
			xlim = range (sf)
		if (missing (ylim) )
			ylim = c (0, 1)
		x = seq (sf)
		y = sf (x)

		if (missing (xlab) )
			xlab = .deflab (sf)
		if (missing (ylab) )
			ylab = "cumprob"

		if (missing (fill.color) )
			fill.color = options$main.fill.color
	}
	else if (is.cqf (sf) )
	{	if (missing (xlim) )
			xlim = c (0, 1)
		if (missing (ylim) )
			ylim = range (sf)
		h0 = ylim [1]
		x = seq (0, 1, length.out=200)
		y = sf (x)

		if (missing (xlab) )
			xlab = "cumprob"
		if (missing (ylab) )
			ylab = .deflab (sf)

		if (missing (fill.color) )
			fill.color = options$main.fill.color.2
	}
	else
		stop ("sf not valid")
	if (is.cqf (sf) || is.cpdc (sf) )
		data = FALSE
	vo = 0.025 * diff (ylim)
	if (data)
		ylim [1] = ylim [2] / (-4.5)
	ylim [2] = ylim [2] + vo

	if (! add)
	{	plot.new ()
		plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
		title (main=main, xlab=xlab, ylab=ylab)
	}

	.pcpd0 (x, y, line.width, line.color, fill.color, h0)
	if (! add)
	{	if (data)
		{	abline (h=0)
			if (sf %$% "is.weighted") .wpoints (sf, ylim [1] + vo, -vo)
			else .xpoints (sf, ylim [1] + vo, -vo, options$main.fill.color.2)
		}
		box ()
		if (axes [1])
			axis (1)
		if (axes [2])
		{	if (is.cqf (sf) )
				axis (2)
			else
			{	yat = pretty (c (0, y), 3)
				yat = yat [yat >= 0]
				axis (2, yat)
			}
		}
	}
}

.xpoints = function (sf, ya, yb, fill.color)
{	x = sf %$% "x"
	y = runif (sf %$% "n", ya, yb)
	points (x, y, pch=16, col=fill.color)
	points (x, y)
}

.wpoints = function (sf, ya, yb)
{	w = sf %$% "w"
	x = sf %$% "x"
	y = runif (sf %$% "n", ya, yb)
	I = order (w)
	.wpoints.2 (w, x, y)
}

.wpoints.2 = function (w, x, y)
{	I = order (w)
	x = x [I]
	y = y [I]
	b = 1 - w [I]
	points (x, y, pch=16, col = rgb (b, b, b) )
	points (x, y)
}
