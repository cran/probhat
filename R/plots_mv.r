#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

setClass ("BVMatrix",
	slots = list (fv="matrix", x="numeric", y="numeric") )

.plot_cfield.ext = function (x, y, fv, ..., add=FALSE, ncontours, contours=TRUE, fb, contour.labels, heatmap=TRUE, xyrel="m",
	.is.cdf, .has.points, .points, .point.color)
{	if (missing (ncontours) )
	{	if (.is.cdf) ncontours = 4
		else ncontours = 6
	}
	if (missing (contour.labels) )
	{	if (.is.cdf) contour.labels = TRUE
		else contour.labels = FALSE
	}
	if (.is.cdf && missing (fb) )
	{	if (ncontours == 1) fb = 0.5
		else if (ncontours == 2) fb = c (0.33, 0.67)
		else if (ncontours == 3) fb = c (0.25, 0.5, 0.75)
		else if (ncontours == 4) fb = c (0.2, 0.4, 0.6, 0.8)
		else
		{	fb = seq (0, 1, length.out = ncontours + 2)
			fb = fb [2:(ncontours + 1)]
			fb = round (fb, 2)
		}
	}

	plot_cfield (x, y, fv, xyrel=xyrel, add=add, contours=FALSE, heatmap=heatmap, ...)
	if (.has.points)
		points (.points [,1], .points [,2], pch=16, col=.point.color)
	plot_cfield (x, y, fv, fb=fb, ncontours=ncontours, contours=contours, contour.labels=contour.labels, heatmap=FALSE, add=TRUE, ...)
}

.plot_surface = function (..., zlab, ref.arrows=FALSE, z.axis, .sf)
{	is.cdf = is.ccdf (.sf)

	if (missing (z.axis) )
	{	if (is.cdf) z.axis=TRUE
		else z.axis=FALSE
	}
	if (missing (zlab) )
	{	if (z.axis)
		{	if (is.cdf) zlab = "Fh"
			else zlab = "fh"
			xlab = .mvvarname (.sf, 0, 1, TRUE, req=TRUE)
			ylab = .mvvarname (.sf, 0, 2, TRUE, req=TRUE)
			if (is.cpdmvc (.sf) )
				constr = " | ..."
			else
				constr = ""
			zlab = paste0 (zlab, " (", xlab, ", ", ylab, constr, ")")
		}
		else
			zlab=""
	}
	plot_surface (..., zlab=zlab, ref.arrows=ref.arrows, z.axis=z.axis)
}

.plotf_isosurface = function (..., ref.arrows=FALSE, maximal=TRUE)
	plotf_isosurface (..., ref.arrows=ref.arrows, maximal=maximal)
.plotf_cfield_3d.ext = function (..., emph="h")
	plotf_cfield_3d (..., emph=emph)

.mvvarname = function (sf, xlab, I, compact=FALSE, ..., req=FALSE)
{	if (req || missing (xlab) )
	{	xlab = (names (sf) ) [I]
		if (is.cpdmvc (sf) && ! compact)
			xlab = paste (xlab, "| ...")
	}
	xlab
}

plot_cpd_bv = function (sf, in3d=FALSE, data, ..., n=30,
	main, xlab, ylab, zlab, xlim, ylim, zlim,
	add=FALSE, point.color)
{	. = attributes (sf)
	if (missing (main) )
		main = ""
	xlab = .mvvarname (sf, xlab, 1)
	ylab = .mvvarname (sf, ylab, 2)

	is.cdf = is.ccdf (sf)
	v = ph4.BVMatrix (sf, xlim, ylim, n=n)
	if (in3d)
	{	if (missing (zlim) )
		{	if (is.cdf) zlim = c (0, 1)
			else zlim = c (0, max (v@fv) )
		}
		.plot_surface (v@x, v@y, v@fv, main=main, xlab=xlab, ylab=ylab, zlab=zlab, zlim=zlim, add=add, ..., .sf=sf)
	}
	else
	{	if (is.cpdmvc (sf) )
			data = FALSE
		else if (missing (data) )
		{	if (add)
				data = FALSE
			else if ( (sf %$% "n") <= 2000)
				data = TRUE
			else
				data = FALSE
		}

		if (data)
		{	ps = (sf %$% "data")$x
			if (missing (point.color) )
			{	options = getOption ("probhat")
				point.color = options$semi.fill.color
			}
		}
		else
			ps = NULL
		
		.plot_cfield.ext (v@x, v@y, v@fv, main=main, xlab=xlab, ylab=ylab, add=add, ...,
			.is.cdf=is.cdf,
			.has.points=data, .points=ps, .point.color=point.color)
	}
}

plot_cpd_tv = function (sf, iso=FALSE, ...,
	main, xlab, ylab, zlab,
	xlim, ylim, zlim,
	z.reverse=FALSE)
{	. = attributes (sf)
	g = function (x, y, z)
		sf (cbind (x, y, z) )
	lims = range (sf)
	if (missing (xlim) ) xlim = lims [1,]
	if (missing (ylim) ) ylim = lims [2,]
	if (missing (zlim) )
	{	zlim = lims [3,]
		if (z.reverse)
			zlim = rev (zlim)
	}
	if (missing (main) ) main = ""
	xlab = .mvvarname (sf, xlab, 1)
	ylab = .mvvarname (sf, ylab, 2)
	zlab = .mvvarname (sf, zlab, 3)

	if (iso)
		.plotf_isosurface (g, xlim, ylim, zlim, ..., main=main, xlab=xlab, ylab=ylab, zlab=zlab)
	else
		.plotf_cfield_3d.ext (g, xlim, ylim, zlim, ..., main=main, xlab=xlab, ylab=ylab, zlab=zlab)
}

ph4.BVMatrix = function (sf, xlim, ylim, ..., n=10)
{	n = rep_len (n, 2)
	lims = range (sf)
	if (missing (xlim) ) xlim = lims [1,]
	if (missing (ylim) ) ylim = lims [2,]
	x = seq (xlim [1], xlim [2], length.out = n [1])
	y = seq (ylim [1], ylim [2], length.out = n [2])
	fv = outer (x, y, .outer.cbind.ext, sf)
	new ("BVMatrix", fv=fv, x=x, y=y)
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
