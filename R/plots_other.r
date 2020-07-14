#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

list.ckernels = function ()
{	ks = vector ("list", 6)
	ks [[1]] = biweight.ckernel ()
	ks [[2]] = truncnorm.ckernel ()
	ks [[3]] = epanechnikov.ckernel ()
	ks [[4]] = triweight.ckernel ()
	ks [[5]] = tricube.ckernel ()
	ks [[6]] = bell.spline ()
	ks
}

.empty.kernel.plot = function (main="", y)
{	plot.new ()
	plot.window (xlim = c (-1.05, 1.05), ylim = c (0, 1.25), xaxs="i", yaxs="i")
	title (main)
	abline (h=y, lty=3)
}

kernel.array = function (ks = list.ckernels (), ..., ref.line=TRUE, colors)
{	n = length (ks)
	if (missing (colors) )
	{	options = getOption ("probhat")
		colors = hcl.colors (n, options$palette, 0.55)
	}
	p0 = par (mfrow = c (n, n), oma = c (2.5, 2.5, 0.25, 0.25), mar = c (0.5, 0.5, 3, 0.5) )
	for (i in 1:n)
	{	for (j in 1:n)
		{	if (ref.line)
				y = ks [[i]]$pdf (0)
			else
				y = NA
			if (i == j)
			{	.empty.kernel.plot (ks [[i]] %$% "name", y)
				plot (ks [[i]], add=TRUE, fill.color = colors [i])
				box ()
				axis (1, at = c (-1, 0, 1) )
				axis (2, at = c (0, 0.5, 1) )
			}
			else if (j > i)
			{	.empty.kernel.plot (,y)
				plot (ks [[i]], add=TRUE, line.color=NA, fill.color = colors [i])
				plot (ks [[j]], add=TRUE, fill.color = colors [j])
				plot (ks [[i]], add=TRUE, fill.color=NA)
				box ()
			}
			else
				plot.new ()
		}
	}
	par (p0)
}

.plot.distribution.set = function (fs, legend, colors, ...)
{	if (is.pmfuv (fs [[1]]) || is.pdfuv (fs [[1]]) || is.ccdfuv (fs [[1]]) )
		.plot.distribution.set.overlay (fs, legend, colors, ...)
	else
		.plot.distribution.set.stacked (fs, colors, ...)
}

.plot.distribution.set.overlay = function (fs, legend, colors, ...)
{	n = length (fs)
	if (missing (colors) )
	{	options = getOption ("probhat")
		colors = hcl.colors (n, options$palette, 0.55)
	}
	xlim = matrix (0, n, 2)
	ymax = numeric (n)
	for (i in 1:n)
	{	xlim [i,] = fs [[i]] %$% "xlim"
		x = seq (fs [[i]], n=200)
		y = fs [[i]](x)
		ymax [i] = max (y)
	}
	xlim = c (min (xlim [,1]), max (xlim [,2]) )
	ylim = c (0, 1.025 * max (ymax) )
	plot (fs [[1]], line.color=NA, fill.color = colors [1], xlim=xlim, ylim=ylim, ...)
	if (n > 1)
	{	for (i in 2:n)
			plot (fs [[i]], line.color=NA, fill.color = colors [i], ..., add=TRUE)
	}
	for (i in 1:n)
		lines (fs [[i]])
	box ()
	if (legend)
	{	legend ("topright",, fs %$% "levnames", colors, bty="n")
	}
}

.plot.distribution.set.stacked = function (fs, nr, nc, colors, ...)
{	n = length (fs)
	if (missing (nr) || missing (nc) )
		p0 = par (mfrow = c (n, 1) )
	else
		p0 = par (mfrow = c (nr, nc) )
	if (missing (colors) )
	{	for (i in 1:n)
			plot (fs [[i]], ...)
	}
	else
	{	for (i in 1:n)
			plot (fs [[i]], area.color = colors [i], ...)
	}
	par (p0)
}

plot.ph3.gset = function (x, ..., legend=TRUE, colors)
	.plot.distribution.set (x, legend, colors, ...)

plot.ph3.mset = function (x, ..., nr, nc, colors)
	.plot.distribution.set.stacked (x, nr, nc, colors, ...)
