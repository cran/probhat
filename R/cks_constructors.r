#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

pdfuv.cks = function (x, spline=TRUE, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
	.cksuv (.pdfuv.cks.eval, "pdfuv", "pdfuv.cks", spline, nc, kernel, bw, smoothness, x, w)

.pdfuv.cks.eval = function (x)
{	. = THAT ()
	if (.$is.spline)
		.$spline.function (x)
	else
		.iterate.uv (.pdfuv.cks.eval.scalar, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
}

cdfuv.cks = function (x, spline=TRUE, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
	.cksuv (.cdfuv.cks.eval, "cdfuv", "cdfuv.cks", spline, nc, kernel, bw, smoothness, x, w)

.cdfuv.cks.eval = function (x, spline=TRUE, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
{	. = THAT ()
	if (.$is.spline)
		.$spline.function (x)
	else
		.iterate.uv (.cdfuv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
}

qfuv.cks = function (x, spline.only, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
{	F = function (x)
	{	. = THAT ()
		.iterate.uv (.cdfuv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
	}
	F.inv = function (y)
	{	. = THAT ()
		.test.y.ok (y)
		.$spline.function (y)
	}
  	.qfuv.cks (F, F.inv, nc, kernel, bw, smoothness, x, w)
}

pdfmv.cks = function (x, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
{	f = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$m, x)
		.iterate.mv (.pdfmv.cks.eval.scalar, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$m, .$x, .$w, u=x)
	}
	.cksmv.2 (f, "pdfmv", "pdfmv.cks", kernel, bw, smoothness, x, w)
}

cdfmv.cks = function (x, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
{	F = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$m, x)
		.iterate.mv (.cdfmv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$m, .$x, .$w, u=x)
	}
	.cksmv.2 (F, "cdfmv", "cdfmv.cks", kernel, bw, smoothness, x, w)
}

pdfc.cks = function (x, spline=TRUE, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
	conditions, preserve.range=FALSE)
{	f = function (x)
	{	. = THAT ()
		if (.$is.spline)
			.$spline.function (x)
		else
			.iterate.uv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksc.2 (f, "pdfuv", "pdfc", "pdfc.cks", spline, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
}

cdfc.cks = function (x, spline=TRUE, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
	conditions, preserve.range=FALSE)
{	F = function (x)
	{	. = THAT ()
		if (.$is.spline)
			.$spline.function (x)
		else
			.iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksc.2 (F, "cdfuv", "cdfc", "cdfc.cks", spline, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
}

qfc.cks = function (x, spline.only, nc=30, kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
	conditions, preserve.range=FALSE)
{	F = function (x)
	{	. = THAT ()
		.iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
	}
	F.inv = function (y)
	{	. = THAT ()
		.test.y.ok (y)
		.$spline.function (y)
	}
	.qfc.cks (F, F.inv, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
}

.cdfc.cks.eval.2 = function (x)
{	. = THAT ()
	.iterate.uv (.cdfc.cks.eval.scalar.2, .$ncon, .$is.weighted, .$conditions, .$kernel$pdf, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
}

.qfc.cks.eval.2 = function (y)
{	. = THAT ()
	.$spline.function (y)
}

pdfmvc.cks = function (x, kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
	conditions, preserve.range=FALSE)
{	f = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$M, x)
		.iterate.mv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksmvc.2 (f, "pdfmv", "pdfmvc", "pdfmvc.cks", preserve.range, conditions, kernel, bw, smoothness, x, w)
}

cdfmvc.cks = function (x, kernel=biweight.kernel, bw, smoothness=0.66, w=NA,
	conditions, preserve.range=FALSE)
{	F = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$M, x)
		.iterate.mv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksmvc.2 (F, "cdfmv", "cdfmvc", "cdfmvc.cks", preserve.range, conditions, kernel, bw, smoothness, x, w)
}

chqf.cks = function (x, spline.only, nc=16, kernel=biweight.kernel, bw, smoothness=0.66, w=NA)
{	chqf.f = function (y)
	{	this.f = THIS ()
		y = .val.u.mv (this.f %$% "m", y)
		x = .chqf.cks.eval (this.f, y)
		colnames (x) = this.f %$% "variable.names"
		x
	}
	.chqf.cks (chqf.f, nc, kernel, bw, smoothness, x, w)
}
