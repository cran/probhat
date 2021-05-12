#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

pdfuv.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	spline=TRUE, bw.method="ph.default", nc=30,
	trtype="local",
	Xlim = cbind (a, b), a=-Inf, b=Inf)
{	.arg.error (...)
	.cksuv (.pdfuv.cks.eval, .CV.pdfuv.cks,
		FALSE, spline, nc, Xlim, bw.method, kernel, trtype, bw, smoothness, x, w)
}

cdfuv.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	spline=TRUE, bw.method="ph.default", nc=30, tail="lower",
	trtype="local",
	Xlim = cbind (a, b), a=-Inf, b=Inf)
{	.arg.error (...)
	.cksuv (.cdfuv.cks.eval, .CV.cdfuv.cks,
		TRUE, spline, nc, Xlim, bw.method, kernel, trtype, bw, smoothness, x, w, tail)
}

qfuv.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", nc=30,
	trtype="local",
	Xlim = cbind (a, b), a=-Inf, b=Inf)
{	.arg.error (...)
	.qfuv.cks (.cdfuv.cks.eval, .qfuv.cks.eval,
		nc, Xlim, bw.method, kernel, trtype, bw, smoothness, x, w, tail="lower")
}

pdfmv.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default",
	Xlim = cbind (a, b), a=-Inf, b=Inf)
{	.arg.error (...)
	.cksmv (.pdfmv.cks.eval, .CV.pdfmv.cks,
		Xlim, bw.method, kernel, bw, smoothness, x, w)
}

cdfmv.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", tail="lower",
	Xlim = cbind (a, b), a=-Inf, b=Inf)
{	.arg.error (...)
	.cksmv (.cdfmv.cks.eval, .CV.cdfmv.cks,
		Xlim, bw.method, kernel, bw, smoothness, x, w, "mv", TRUE, tail=tail)
}

pdfc.cks = function (x, ..., conditions, w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	spline=TRUE, bw.method="ph.default", nc=30,
	Xlim = cbind (a, b), a=-Inf, b=Inf,
	preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE, warning=TRUE)
{	.arg.error (...)
	.cksc.2 (.pdfc.cks.eval, .CV.pdfc.cks,
		FALSE, spline, nc, preserve.range, conditions, Xlim, bw.method, kernel, bw, smoothness, x, w,
		as.cset, as.list.cset, warning)
}

cdfc.cks = function (x, ..., conditions, w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	spline=TRUE, bw.method="ph.default", nc=30, tail="lower",
	Xlim = cbind (a, b), a=-Inf, b=Inf,
	preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE, warning=TRUE)
{	.arg.error (...)
	.cksc.2 (.cdfc.cks.eval, .CV.cdfc.cks,
		TRUE, spline, nc, preserve.range, conditions, Xlim, bw.method, kernel, bw, smoothness, x, w,
		as.cset, as.list.cset, warning, tail)
}

qfc.cks = function (x, ..., conditions, w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", nc=30,
	Xlim = cbind (a, b), a=-Inf, b=Inf,
	preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE, warning=TRUE)
{	.arg.error (...)
	.qfc.cks (.cdfc.cks.eval, .qfc.cks.eval,
		nc, preserve.range, conditions, Xlim, bw.method, kernel, bw, smoothness, x, w,
		as.cset, as.list.cset, warning, tail="lower")
}

pdfmvc.cks = function (x, ..., conditions, w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default",
	Xlim = cbind (a, b), a=-Inf, b=Inf,
	preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE, warning=TRUE)
{	.arg.error (...)
	.cksmvc.2 (.pdfmvc.cks.eval, .CV.pdfmvc.cks, FALSE,
		preserve.range, conditions, Xlim, bw.method, kernel, bw, smoothness, x, w,
		as.cset, as.list.cset, warning)
}

cdfmvc.cks = function (x, ..., conditions, w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", tail="lower",
	Xlim = cbind (a, b), a=-Inf, b=Inf,
	preserve.range=FALSE, as.cset=FALSE, as.list.cset=FALSE, warning=TRUE)
{	.arg.error (...)
	.cksmvc.2 (.cdfmvc.cks.eval, .CV.cdfmvc.cks, TRUE,
		preserve.range, conditions, Xlim, bw.method, kernel, bw, smoothness, x, w,
		as.cset, as.list.cset, warning, tail)
}

chqf.cks = function (x, ..., w,
	bw, smoothness=1,
	kernel=BIWEIGHT.CKERNEL,
	bw.method="ph.default", nc=16)
{	.arg.error (...)
	.chqf.cks (.chqf.cks.eval,
		nc, bw.method, kernel, bw, smoothness, x, w)
}
