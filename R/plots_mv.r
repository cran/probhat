#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.plot_cfield.ext = function (x, y, fv, ..., add=FALSE,contours=TRUE, heatmap=TRUE, xyrel="m", .points)
{	options = options = getOption ("probhat")
	plot_cfield (x, y, fv, xyrel=xyrel, add=add, contours=FALSE, heatmap=heatmap, ...)
	if (! is.null (.points) )
		points (.points [,1], .points [,2], pch=16, col=options$semi.fill.color)
	plot_cfield (x, y, fv, contours=contours, heatmap=FALSE, add=TRUE, ...)
}

.plotf_cfield_3d.ext = function (..., emph="h") plotf_cfield_3d (..., emph=emph)

.mvvarname = function (sf, xlab, I)
{	if (missing (xlab) )
	{	xlab = (sf %$% "variable.names") [I]
		if (is.cpdmvc (sf) )
			xlab = paste ("conditional", xlab)
	}
	xlab
}

plot_cpd_bv = function (sf, in3d=FALSE, data=FALSE, ..., all=FALSE, n=30,
	main, xlab, ylab, xlim, ylim, zlim)
{	. = attributes (sf)
	Jx = sf %$% "m" - 1
	Jy = Jx + 1
	if (missing (main) )
		main = ""
	xlab = .mvvarname (sf, xlab, Jx)
	ylab = .mvvarname (sf, ylab, Jy)
	if (is.cpdmvc (sf) )
	{	all = FALSE
		data = FALSE
	}
	if (all)
	{	f1 = .uv.from.mv.cks (sf, Jx)
		f2 = .uv.from.mv.cks (sf, Jy)
		mar = c (4.75, 4.75, 0.75, 0.75)
		p0 = par (cex=0.65, mar=mar, mfcol = c (2, 2) )
		plot (f1, data)
		plot (f2, data)
		plot_cpd_bv (sf, FALSE, data, ..., n=n)
		plot_cpd_bv (sf, TRUE, ..., n=n)
		par (p0)
	}
	else
	{	if (missing (xlim) )
			xlim = .$xlim [1,]
		if (missing (ylim) )
			ylim = .$xlim [2,]
		n = rep_len (n, 2)
		x = seq (xlim [1], xlim [2], length.out = n [1])
		y = seq (ylim [1], ylim [2], length.out = n [2])
		z = outer (x, y, .outer.cbind.ext, sf)
		if (in3d)
		{	if (missing (zlim) )
			{	if (is.ccdf (sf) )
					zlim = c (0, 1)
				else
					zlim = c (0, max (z) )
			}
			plot_surface (x, y, z, main=main, xlab=xlab, ylab=ylab, zlim=zlim, ...)

		}
		else
		{	if (data)	p = sf %$% "x"
			else p = NULL
			.plot_cfield.ext (x, y, z, main=main, xlab=xlab, ylab=ylab, ..., .points=p)
		}
	}
}

plot_cpd_tv = function (sf, ...,
	main, xlab, ylab, zlab,
	xlim, ylim, zlim,
	reverse.z=FALSE)
{	. = attributes (sf)
	g = function (x, y, z)
		sf (cbind (x, y, z) )
	Jx = sf %$% "m" - 2
	Jy = Jx + 1
	Jz = Jx + 2
	if (missing (xlim) ) xlim = .$xlim [1,]
	if (missing (ylim) ) ylim = .$xlim [2,]
	if (missing (zlim) )
	{	zlim = .$xlim [3,]
		if (reverse.z) zlim = rev (zlim)
	}
	if (missing (main) ) main = ""
	xlab = .mvvarname (sf, xlab, Jx)
	ylab = .mvvarname (sf, ylab, Jy)
	zlab = .mvvarname (sf, zlab, Jz)
	
	.plotf_cfield_3d.ext (g, xlim, ylim, zlim, ..., main=main, xlab=xlab, ylab=ylab, zlab=zlab)
}

.outer.cbind.ext = function (x, y, f)
	f (cbind (x, y) )

.outer.3 = function (f, n, x, y, z)
{	x2 = rep (x, times = n [2] * n [3])
	y2 = rep (y, each = n [1], times = n [3])
	z2 = rep (z, each = n [1] * n [2])
	fv = f (cbind (x2, y2, z2) )
	array (fv, n)
}
