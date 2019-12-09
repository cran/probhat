#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.test.bv = function (M)
{	if (M != 2)
		stop ("can only plot multivariate models, if bivariate")
}

.plot.cksmv = function (f, use.plot3d=FALSE,
	main, xlab, ylab,
	npoints=30, xlim, ylim, ...,
	all=FALSE, contrast=-0.8, cex=0.65)
{	. = attributes (f)
	Jx = f %$% "m" - 1
	Jy = Jx + 1
	if (missing (main) )
		main = ""
	if (missing (xlab) )
		xlab = (f %$% "variable.names") [Jx]
	if (missing (ylab) )
		ylab = (f %$% "variable.names") [Jy]
	if (all)
	{	f1 = .uv.from.mv.cks (f, Jx)
		f2 = .uv.from.mv.cks (f, Jy)
		mar = c (4.75, 4.75, 0.75, 0.75)
		p0 = par (cex=cex, mar=mar)
		layout (matrix (c (1, 3, 4, 2, 5, 5), 3, 2),, c (1, 0.5, 0.5) )
		.plot.cksmv (f, FALSE, main, xlab, ylab, contrast=contrast)
		.plot.cksmv (f, TRUE, main, xlab, ylab, npoints)
		p1 = par (mar = c (2.75, 2.75, 0.75, 0.75), mgp = c (0.5, 0, 0) )
		plot (f1, axis.ticks=FALSE)
		plot (f2, axis.ticks=FALSE)
		par (p1)
		.plot.points (f %$% "is.weighted", (f %$% "x") [,Jx], (f %$% "x") [,Jy], f %$% "w", xlab, ylab)
		par (p0)
	}
	else
	{	if (missing (xlim) )
			xlim = .$xlim [1,]
		if (missing (ylim) )
			ylim = .$xlim [2,]
		x = seq (xlim [1], xlim [2], length.out=npoints)
		y = seq (ylim [1], ylim [2], length.out=npoints)
		z = outer (x, y, .outer.cbind.ext, f)
		if (use.plot3d)
			plot3d.surface (,,z, main=main, xlab=xlab, ylab=ylab, ...)
		else
			plot2d.contour (x, y, z, main=main, xlab=xlab, ylab=ylab, contrast=contrast, ...)
	}
}

plot.pdfmv.cks = function (x, use.plot3d=FALSE, main, xlab, ylab, npoints=30, ..., all=FALSE)
{	.test.bv (x %$% "m")
	.plot.cksmv (x, use.plot3d, main, xlab, ylab, npoints, ..., all=all)
}

plot.cdfmv.cks = function (x, use.plot3d=FALSE, main, xlab, ylab, npoints=30, ..., all=FALSE)
{	.test.bv (x %$% "m")
	.plot.cksmv (x, use.plot3d, main, xlab, ylab, npoints, ..., all=all)
}

plot.pdfmvc.cks = function (x, use.plot3d=FALSE, main, xlab, ylab, npoints=30, ...)
{	.test.bv (x %$% "M")
	.plot.cksmv (x, use.plot3d, main, xlab, ylab, npoints, ..., all=FALSE)
}

plot.cdfmvc.cks = function (x, use.plot3d=FALSE, main, xlab, ylab, npoints=30, ...)
{	.test.bv (x %$% "M")
	.plot.cksmv (x, use.plot3d, main, xlab, ylab, npoints, ..., all=FALSE)
}

plot.chqf = function (x, ...)
	stop ("can't plot chained quantile function")

.outer.cbind.ext = function (x, y, f)
	f (cbind (x, y) )
