#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.pdfuv.cks.eval = function (x)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc, .$XLIM)
	if (.$is.spline)
		.$spline.function (x)
	else
	{	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
		y = .iterate.uv (.pdfuv.cks.eval.scalar, .$is.weighted, .$kernel@f, .$bw, data$n, data$x, data$w, u=x)
		.scale.val (y, .$.any.trunc, .$.scalef)
	}
}

.cdfuv.cks.eval = function (x, ..., .ignore.trunc=FALSE)
{	. = .THAT ()
	if (! .ignore.trunc)
		x = .val.u.uv (x, .$.any.trunc, .$.is.trunc, .$XLIM)
	if (.$is.spline)
		.$spline.function (x)
	else
	{	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
		y = .iterate.uv (.cdfuv.cks.eval.scalar, .$is.weighted, .$kernel@F, .$bw, data$n, data$x, data$w, .$.low, u=x)
		if (! .ignore.trunc)
		{	if (.$.any.trunc.lower)
				y = y - .$.const.cdf.lower
			y = .scale.val (y, .$.any.trunc, .$.scalef)
		}
		y
	}
}

.qfuv.cks.eval = function (p)
{	. = .THAT ()
	.test.y.ok (p)
	.$spline.function (p)
}

.pdfmv.cks.eval = function (x)
{	. = .THAT ()
	x = .val.u.mv (.$m, x, .$.any.trunc, .$.is.trunc, .$XLIM)
	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
	y = .iterate.mv (.pdfmv.cks.eval.scalar, .$is.weighted, .$kernel@f, .$m, .$bw, data$n, data$x, data$w, u=x)
	.scale.val (y, .$.any.trunc, .$.scalef)
}

.cdfmv.cks.eval = function (x, ..., .ignore.trunc=FALSE)
{	this = .THIS ()
	. = .THAT ()
	if (! .ignore.trunc)
		x = .val.u.mv (.$m, x, .$.any.trunc, .$.is.trunc, .$XLIM)
	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
	y = .iterate.mv (.cdfmv.cks.eval.scalar, .$is.weighted, .$kernel@F, .$m, .$bw, data$n, data$x, data$w, .$.low, u=x)
	if (! .ignore.trunc)
	{	if (.$.any.trunc.lower)
		{	for (i in seq_len (nrow (x) ) )
				y [i] = y [i] - .cdfv.lower.side (this, .$.low, .$.is.trunc.lower, .$m, .$data$xlim, hi = x [i,])
		}
		y = .scale.val (y, .$.any.trunc, .$.scalef)
	}
	y
}

.pdfc.cks.eval = function (x)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc [.$m,], .$XLIM [.$m,])
	if (.$is.spline)
		.$spline.function (x)
	else
	{	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
		y = .iterate.uv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel@f, .$bw, data$n, data$x, data$w, u=x)
		.scale.val (y, .$.any.trunc, .$.scalef)
	}
}

.cdfc.cks.eval = function (x, ..., .ignore.trunc=FALSE)
{	. = .THAT ()
	if (! .ignore.trunc)
		x = .val.u.uv (x, .$.any.trunc, .$.is.trunc [.$m,], .$XLIM [.$m,])
	if (.$is.spline)
		.$spline.function (x)
	else
	{	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
		y = .iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel@F, .$bw, data$n, data$x, data$w, .$.low, u=x)
		if (! .ignore.trunc)
		{	if (.$.any.trunc.lower)
				y = y - .$.const.cdf.lower
			y = .scale.val (y, .$.any.trunc, .$.scalef)
		}
		y
	}
}

.qfc.cks.eval = function (p)
{	. = .THAT ()
	.test.y.ok (p)
	.$spline.function (p)
}

.pdfmvc.cks.eval = function (x)
{	. = .THAT ()
	J = (.$ncon + 1):(.$m)
	x = .val.u.mv (.$M, x, .$.any.trunc, .$.is.trunc [J,, drop=FALSE], .$XLIM [J,, drop=FALSE])
	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
	y = .iterate.mv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel@f, .$bw, data$n, data$x, data$w, u=x)
	.scale.val (y, .$.any.trunc, .$.scalef)
}

.cdfmvc.cks.eval = function (x, ..., .ignore.trunc=FALSE)
{	this = .THIS ()
	. = .THAT ()
	J = (.$ncon + 1):(.$m)
	if (! .ignore.trunc)
		x = .val.u.mv (.$M, x, .$.any.trunc, .$.is.trunc [J,, drop=FALSE], .$XLIM [J,, drop=FALSE])
	data = .select.bdata (.$.any.trunc, .$data, .$.xpnd)
	y = .iterate.mv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel@F, .$bw, data$n, data$x, data$w, .$.low, u=x)
	if (! .ignore.trunc)
	{	if (.$.any.trunc.lower)
		{	for (i in seq_len (nrow (x) ) )
				y [i] = y [i] - .cdfv.lower.side (this, .$.low, .$.is.trunc.lower, .$M, .$data$xlim [J,, drop=FALSE],
					hi = x [i,])
	
		}
		y = .scale.val (y, .$.any.trunc, .$.scalef)
	}
	y
}

.chqf.cks.eval = function (p)
{	this.f = .THIS ()
	p = .val.u.mv (this.f %$% "m", p)
	x = .chqf.cks.eval.ext (this.f, p)
	colnames (x) = this.f %$% "xnames"
	x
}

.cdfc4chqf.cks.eval = function (x)
{	. = .THAT ()
	data = .$data
	.iterate.uv (.cdfc4chqf.cks.eval.scalar, .$ncon, .$is.weighted, .$conditions, .$kernel@f, .$kernel@F, .$bw, data$n, data$x, data$w, u=x)
}

.qfc4chqf.cks.eval = function (p)
{	. = .THAT ()
	.$spline.function (p)
}

.scale.val = function (y, trunc, k)
{	if (trunc) k * y
	else y
}

.select.bdata = function (trunc, data, xpnd)
{	if (trunc) xpnd
	else data
}
