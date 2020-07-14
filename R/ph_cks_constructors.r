#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

pdfuv.cks = function (x, ..., spline=TRUE, bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA)
	.cksuv (.pdfuv.cks.eval, .CV.pdfuv.cks, FALSE, spline, nc, bw.method, kernel, bw, smoothness, x, w)

.pdfuv.cks.eval = function (x)
{	. = THAT ()
	if (.$is.spline)
		.$spline.function (x)
	else
		.iterate.uv (.pdfuv.cks.eval.scalar, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
}

cdfuv.cks = function (x, ..., spline=TRUE, bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA)
	.cksuv (.cdfuv.cks.eval, .CV.cdfuv.cks, TRUE, spline, nc, bw.method, kernel, bw, smoothness, x, w)

.cdfuv.cks.eval = function (q)
{	. = THAT ()
	if (.$is.spline)
		.$spline.function (q)
	else
		.iterate.uv (.cdfuv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
}

qfuv.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA)
{	F = function (q)
	{	. = THAT ()
		.iterate.uv (.cdfuv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
	}
	F.inv = function (p)
	{	. = THAT ()
		.test.y.ok (p)
		.$spline.function (p)
	}
  	.qfuv.cks (F, F.inv, nc, bw.method, kernel, bw, smoothness, x, w)
}

pdfmv.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, bw, smoothness=1, w=NA)
{	f = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$m, x)
		.iterate.mv (.pdfmv.cks.eval.scalar, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$m, .$x, .$w, u=x)
	}
	.cksmv (f, .CV.pdfmv.cks, bw.method, kernel, bw, smoothness, x, w)
}

cdfmv.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, bw, smoothness=1, w=NA)
{	F = function (q)
	{	. = THAT ()
		q = .val.u.mv (.$m, q)
		.iterate.mv (.cdfmv.cks.eval.scalar, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$m, .$x, .$w, u=q)
	}
	.cksmv (F, .CV.cdfmv.cks, bw.method, kernel, bw, smoothness, x, w)
}

pdfc.cks = function (x, ..., spline=TRUE, bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA,
	conditions, preserve.range=FALSE, warning=TRUE)
{	f = function (x)
	{	. = THAT ()
		if (.$is.spline)
			.$spline.function (x)
		else
			.iterate.uv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksc.2 (f, .CV.pdfc.cks, FALSE, spline, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
}

cdfc.cks = function (x, ..., spline=TRUE, bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA,
	conditions, preserve.range=FALSE, warning=TRUE)
{	F = function (q)
	{	. = THAT ()
		if (.$is.spline)
			.$spline.function (q)
		else
			.iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
	}
	.cksc.2 (F, .CV.cdfc.cks, TRUE, spline, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
}

qfc.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, nc=30, bw, smoothness=1, w=NA,
	conditions, preserve.range=FALSE, warning=TRUE)
{	F = function (q)
	{	. = THAT ()
		.iterate.uv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
	}
	F.inv = function (p)
	{	. = THAT ()
		.test.y.ok (p)
		.$spline.function (p)
	}
	.qfc.cks (F, F.inv, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
}

.cdfc.cks.eval.2 = function (q)
{	. = THAT ()
	.iterate.uv (.cdfc.cks.eval.scalar.2, .$ncon, .$is.weighted, .$conditions, .$kernel$pdf, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
}

.qfc.cks.eval.2 = function (p)
{	. = THAT ()
	.$spline.function (p)
}

pdfmvc.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, bw, smoothness=1, w=NA,
	conditions, preserve.range=FALSE, warning=TRUE)
{	f = function (x)
	{	. = THAT ()
		x = .val.u.mv (.$M, x)
		.iterate.mv (.pdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$pdf, .$bw, .$n, .$x, .$w, u=x)
	}
	.cksmvc.2 (f, .CV.pdfmvc.cks, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
}

cdfmvc.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, bw, smoothness=1, w=NA,
	conditions, preserve.range=FALSE, warning=TRUE)
{	F = function (q)
	{	. = THAT ()
		q = .val.u.mv (.$M, q)
		.iterate.mv (.cdfc.cks.eval.scalar, .$.constant, .$.v, .$M, .$ncon, .$is.weighted, .$kernel$cdf, .$bw, .$n, .$x, .$w, u=q)
	}
	.cksmvc.2 (F, .CV.cdfmvc.cks, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
}

chqf.cks = function (x, ..., bw.method="ph.default", kernel=biweight.ckernel, nc=16, bw, smoothness=1, w=NA)
{	chqf.f = function (p)
	{	this.f = THIS ()
		p = .val.u.mv (this.f %$% "m", p)
		x = .chqf.cks.eval (this.f, p)
		colnames (x) = this.f %$% "variable.names"
		x
	}
	.chqf.cks (chqf.f, nc, bw.method, kernel, bw, smoothness, x, w)
}
