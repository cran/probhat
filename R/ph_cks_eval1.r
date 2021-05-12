#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

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
	{	data = .select.bdata (.$.any.trunc, .$trtype, .$data, .$.xpnd)
		y = .iterate.uv (.pdfuv.cks.eval.scalar, .$.internal.isw, .$kernel@f, .$bw, data$n, data$x, .$.internalw, u=x)
		.scale.val (y, .$trtype, .$.any.trunc, .$.scalef)
	}
}

.cdfuv.cks.eval = function (x, ...)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc, .$XLIM)
	if (.$is.spline)
		.$spline.function (x)
	else
	{	data = .select.bdata (.$.any.trunc, .$trtype, .$data, .$.xpnd)
		y = .iterate.uv (.cdfuv.cks.eval.scalar, .$.internal.isw, .$kernel@F, .$bw,
			data$n, data$x, .$.internalw, .$.low, .$.constv, u=x)
		if (.$trtype != "local" && .$.any.trunc.lower)
			y = y - .$.const.cdf.lower
		.scale.val (y, .$trtype, .$.any.trunc, .$.scalef)
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
	.iterate.mv (.pdfmv.cks.eval.scalar, .$.internal.isw, .$kernel@f, .$m, .$bw,
		.$data$n, .$data$x, .$.internalw, u=x)
}

.cdfmv.cks.eval = function (x, ...)
{	. = .THAT ()
	x = .val.u.mv (.$m, x, .$.any.trunc, .$.is.trunc, .$XLIM)
	.iterate.mv (.cdfmv.cks.eval.scalar, .$.internal.isw, .$kernel@F, .$m, .$bw,
		.$data$n, .$data$x, .$.internalw, .$.low, .$.constv, u=x)
}

.pdfc.cks.eval = function (x)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc [.$m,], .$XLIM [.$m,])
	if (.$is.spline)
		.$spline.function (x)
	else
	{	.iterate.uv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$.internal.isw, .$kernel@f, .$bw,
			.$data$n, .$data$x, .$.internalw, u=x)
	}
}

.cdfc.cks.eval = function (x)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc [.$m,], .$XLIM [.$m,])
	if (.$is.spline)
		.$spline.function (x)
	else
	{	.iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$.internal.isw, .$kernel@F, .$bw,
			.$data$n, .$data$x, .$.internalw, .$.low, .$.constv, u=x)
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
	.iterate.mv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$.internal.isw, .$kernel@f, .$bw,
		.$data$n, .$data$x, .$.internalw, u=x)
}

.cdfmvc.cks.eval = function (x)
{	. = .THAT ()
	J = (.$ncon + 1):(.$m)
	x = .val.u.mv (.$M, x, .$.any.trunc, .$.is.trunc [J,, drop=FALSE], .$XLIM [J,, drop=FALSE])
	.iterate.mv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$.internal.isw, .$kernel@F, .$bw,
		.$data$n, .$data$x, .$.internalw, .$.low, .$.constv, u=x)
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

.scale.val = function (y, trtype, trunc, k)
{	if (trunc && (trtype != "local") ) k * y
	else y
}

.select.bdata = function (trunc, trtype, data, xpnd)
{	if (trunc && trtype == "reflect") xpnd
	else data
}
