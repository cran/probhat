#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.pdpd0 = function (x, h, combine=FALSE, space=0, line.width=1, line.color="black", fill.color=NA, hmin=0)
{	n = length (h)
	xb = .ph.bars.x (n, x)
	if (combine)
	{	N = 2 * n + 2
		u = v = rep (hmin, N)
		u [1] = xb [1]
		for (i in 1:n)
		{	I = 2 * i
			u [I] = xb [i]
			u [I + 1] = xb [i + 1]
			v [I] = v [I + 1] = h [i]
		}
		u [N] = xb [n + 1]
		polygon (u, v, border=NA, col=fill.color)
		lines (u, v, lwd=line.width, col=line.color)
	}
	else if (space == 0)
	{	for (i in 1:n)
		{	px = xb [c (i, i, i + 1, i + 1)]
			py = c (hmin, h [i], h [i], hmin)
			polygon (px, py, border=NA, col=fill.color)
			lines (px [1:3], py [1:3], lwd=line.width, col=line.color)
			if (i < n && h [i] > h [i + 1])
				lines (px [3:4], c (h [i], h [i + 1]), lwd=line.width, col=line.color)
		}
		px = xb [c (n + 1, n + 1)]
		py = c (hmin, h [n])
		lines (px, py, lwd=line.width, col=line.color)
		lines (range (xb), c (hmin, hmin), lwd=line.width, col=line.color)
	}
	else
	{	space = abs (xinch (space * 0.0394) )
	
		hs = space / 2
		for (i in 1:n)
		{	px = xb [c (i, i, i + 1, i + 1)]
			px [1:2] = px [1:2] + hs
			px [3:4] = px [3:4] - hs
			py = c (hmin, h [i], h [i], hmin)
			polygon (px, py, lwd=line.width, border=line.color, col=fill.color)
		}
	}
}

.pdpd1 = function (x, h, line.color="black", fill.color=NA, hmin=0)
{	n = length (x)
	for (i in 1:n)
	{	px = c (x [i] - 0.5, x [i] - 0.5, x [i] + 0.5, x [i] + 0.5)
		py = c (hmin, h [i], h [i], hmin)
		polygon (px, py, border=line.color, col=fill.color)
	}
}

plot_dpd = function (sf, data=FALSE, ...,
	main, xlab, ylab,
	xlim, ylim,
	add=FALSE, axes=TRUE, combine = is.dks (sf), freq=FALSE, n, space=0,
	line.width, line.color, fill.color)
{	axes = rep_len (axes, 2)
	options = getOption ("probhat")
	if (missing (main) )
		main = ""
	if (missing (line.width) )
		line.width = options$main.line.width
	if (missing (line.color) )
		line.color = options$main.line.color

	if (missing (space) )
	{	if (is.pmf (sf) && is.cat (sf) )
			space = 2
		else
			space = 0
	}

	if (is.pmf (sf) )
	{	h0 = 0
		if (missing (xlim) )
			xlim = range (sf) + c (-0.5, 0.5)
		x = xb = seq (sf)
		y = sf (x, freq=freq, n=n)
		if (missing (ylim) )
			ylim = c (0, max (y) )
		y2 = sf %$% ".probs"

		if (missing (xlab) )
			xlab = .deflab (sf)
		if (missing (ylab) )
		{	if (freq) ylab = "frequency"
			else ylab = "mass"
		}

		if (missing (fill.color) )
			fill.color = options$main.fill.color
	}
	else if (is.dcdf (sf) )
	{	h0 = 0
		if (missing (xlim) )
			xlim = range (sf) + c (-0.5, 0.5)
		x = xb = seq (sf)
		y = sf (x, freq=freq, n=n)
		if (missing (ylim) )
			ylim = c (0, max (y) )
		y2 = sf %$% ".PROBS"

		if (missing (xlab) )
			xlab = .deflab (sf)
		if (missing (ylab) )
			ylab = "cumprob"

		if (missing (fill.color) )
			fill.color = options$main.fill.color
	}
	else if (is.dqf (sf) )
	{	if (missing (xlim) )
			xlim = c (0, 1)
		if (missing (ylim) )
			ylim = range (sf) + c (-1, 0)
		h0 = ylim [1]
		x = seq (sf, TRUE)
		xb = seq (sf, TRUE, midpoints=FALSE)
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
	if (is.dqf (sf) || is.dpdc (sf) )
		data = FALSE
	vo = 0.025 * diff (ylim)
	if (data)
		ylim [1] = ylim [2] / (-4.5)
	ylim [2] = ylim [2] + vo

	if (! add)
	{	plot.new ()
		if (is.pmf (sf) && is.cat (sf) )
			plot.window (xlim=xlim, ylim=ylim, yaxs="i")
		else
			plot.window (xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
		title (main=main, xlab=xlab, ylab=ylab)
	}
	.pdpd0 (xb, y, combine, space, line.width, line.color, fill.color, h0)
	if (! add)
	{	if (data)
		{	x2 = sf %$% "x"
			h2 = sf %$% "h"
			y2 = ylim [1] + (h2) * -0.9 * ylim [1] / max (h2)

			abline (h=0)
			.pdpd1 (x2, y2, line.color, options$main.fill.color.2, ylim [1])
		}
		box ()
		if (is.dqf (sf) )
		{	if (axes [1]) axis (1)
			if (axes [2])
			{	if (is.cat (sf) )
					axis (2, 1:(sf %$% "nlevels") - 0.075, sf %$% "levels")
				else
					axis (2)
			}
		}
		else
		{	if (is.dks (sf) )
				blabs = .dpd.blabs (x, sf %$% "x")
			else
				blabs = sf %$% "levels"
			if (axes [1]) axis (1, x, blabs)
			if (axes [2])
			{	yat = pretty (c (0, y), 3)
				yat = yat [yat >= 0]
				axis (2, yat)
			}
		}
	}
}

.ph.bars.x = function (n, x)
{	if (missing (x) )
		xb = c (0.5, 1:n + 0.5)
	else
	{	nx = length (x)
		if (nx == n)
		{	dx = mean (diff (x) )
			xb = (x [-1] + x [-n]) / 2
			xb = c (xb [1] - dx, xb, xb [n - 1] + dx)
		}
		else if (nx == n + 1)
			xb = x
		else
			stop ("length (x) not n or n + 1")
	}
	xb
}

.dpd.blabs = function (u, x)
{	rnames = names (x)
	if (is.null (rnames) )
		u
	else
	{	blabs = rep ("", length (u) )
		blabs [match (x, u)] = rnames
		blabs
	}
}
