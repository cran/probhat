#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

plot_2x2 = function (f1, f2, f3, f4, main.1="", main.2="", main.3="", main.4="", ...)
{	p0 = par (mfrow = c (2, 2), mar = c (2.75, 2.75, 3.25, 0.75), cex=0.65)
	.plot_2x2.ext (f1, main.1, ...)
	.plot_2x2.ext (f2, main.2, ...)
	.plot_2x2.ext (f3, main.3, ...)
	.plot_2x2.ext (f4, main.4, ...)
	par (p0)
}

.plot_2x2.ext = function (f, main, ...)
{	if (missing (f) )
		plot.new ()
	else
		plot (f, ..., main=main, xlab="", ylab="")
}

ckernels = function ()
	c ("biweight.kernel", "truncnorm.kernel", "epanechnikov.kernel",
		"triweight.kernel", "tricube.kernel", "bell.spline")

kernel.array = function (which = ckernels (), reference.line=TRUE, colors)
{	n = length (which)
	ks = vector ("list", n)
	strs = paste (which, "()")
	for (i in 1:n)
		ks [[i]] = eval (parse (text = strs [i]) )
	if (missing (colors) )
	{	options = getOption ("probhat")
		colors = hcl.colors (n, options$palette, 0.55)
	}
	p0 = par (mfrow = c (n, n), oma = c (2.5, 2.5, 0.25, 0.25), mar = c (0.5, 0.5, 3, 0.5) )
	for (i in 1:n)
	{	for (j in 1:n)
		{	if (reference.line)
				y = ks [[i]]$pdf (0)
			else
				y = NA
			if (i == j)
			{	.empty.kernel.plot (ks [[i]] %$% "name", y)
				plot (ks [[i]], add=TRUE, area.color = colors [i])
				box ()
				axis (1, at = c (-1, 0, 1) )
				axis (2, at = c (0, 0.5, 1) )
			}
			else if (j > i)
			{	.empty.kernel.plot (,y)
				plot (ks [[i]], add=TRUE, line.color=NA, area.color = colors [i])
				plot (ks [[j]], add=TRUE, area.color = colors [j])
				plot (ks [[i]], add=TRUE, area.color=NA)
				box ()
			}
			else
				plot.new ()
		}
	}
	par (p0)
}

.plot.distribution.set = function (fs, main, xlab, ylab, legend, colors, ...)
{	if (is.pmfuv (fs [[1]]) || is.pdfuv (fs [[1]]) || is.cdfuv (fs [[1]]) )
		.plot.distribution.set.overlay (fs, main, xlab, ylab, legend, colors, ...)
	else
		.plot.distribution.set.stacked (fs, main, xlab, ylab, colors, ...)
}

.plot.distribution.set.overlay = function (fs, main, xlab, ylab, legend, colors, ...)
{	n = length (fs)
	if (missing (colors) )
	{	options = getOption ("probhat")
		colors = hcl.colors (n, options$palette, 0.55)
	}
	xlim = matrix (0, n, 2)
	ymax = numeric (n)
	for (i in 1:n)
	{	xlim [i,] = fs [[i]] %$% "xlim"
		x = seq (fs [[i]], 200)
		y = fs [[i]](x)
		ymax [i] = max (y)
	}
	xlim = c (min (xlim [,1]), max (xlim [,2]) )
	ylim = c (0, 1.025 * max (ymax) )
	plot (fs [[1]], xlab = fs %$% "varname", line.color=NA, area.color = colors [1], xlim=xlim, ylim=ylim, ...)
	if (n > 1)
	{	for (i in 2:n)
			plot (fs [[i]], line.color=NA, area.color = colors [i], ..., add=TRUE)
	}
	for (i in 1:n)
		lines (fs [[i]])
	box ()
	if (legend)
	{	legend ("topright",, fs %$% "levnames", colors, bty="n")
	}
}

.plot.distribution.set.stacked = function (fs, main, xlab, ylab, colors, ...)
{	n = length (fs)
	p0 = par (mfrow = c (n, 1) )
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

plot.marginal.set = function (x, main, xlab, ylab, colors, ...)
	.plot.distribution.set.stacked (x, main, xlab, ylab, colors, ...)

plot.categorical.set = function (x, main, xlab, ylab, legend=TRUE, colors, ...)
	.plot.distribution.set (x, main, xlab, ylab, legend, colors, ...)

plot.conditional.set = function (x, main, xlab, ylab, legend=TRUE, colors, ...)
	.plot.distribution.set (x, main, xlab, ylab, legend, colors, ...)
