#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.k = function (kpdf, bw, x, u)
{	dist = u - x
	2 / bw * kpdf (2 / bw * dist)
}

.K = function (kcdf, bw, x, u, low=TRUE, constv=NA)
{	dist = u - x
	p = kcdf (2 / bw * dist)
	if (! low)
		p = 1 - p
	if (! is.na (constv [1]) )
		p = p - constv
	p
}

.sumk = function (is.weighted, n, w, y)
{	if (is.weighted)
		sum (w * y)
	else
		sum (y) / n
}

.pdfuv.cks.eval.scalar = function (is.weighted, kpdf, bw, n, x, w, u)
{	y = .k (kpdf, bw, x, u)
	.sumk (is.weighted, n, w, y)
}

.cdfuv.cks.eval.scalar = function (is.weighted, kcdf, bw, n, x, w, low, constv, u)
{	y = .K (kcdf, bw, x, u, low, constv)
	.sumk (is.weighted, n, w, y)
}

.pdfmv.cks.eval.scalar = function (is.weighted, kpdf, m, bw, n, x, w, u)
{	y = rep (1, n)
	for (j in 1:m)
		y = y * .k (kpdf, bw [j], x [,j], u [j])
	.sumk (is.weighted, n, w, y)
}

.cdfmv.cks.eval.scalar = function (is.weighted, kcdf, m, bw, n, x, w, low, constv, u)
{	y = rep (1, n)
	for (j in 1:m)
		y = y * .K (kcdf, bw [j], x [,j], u [j], low [j], constv [,j])
	.sumk (is.weighted, n, w, y)
}

#used in both c/mvc
.precompute.cksc.v = function (ncon, conditions, kpdf, bw, n, x)
{	v = rep (1, n)
	for (j in seq_len (ncon) )
		v = v * .k (kpdf, bw [j], x [,j], conditions [j])
	v
}

#used in both c/mvc
.pdfc.cks.eval.scalar = function (constant, v, M, ncon, is.weighted, kpdf, bw, n, x, w, u)
{	for (j in 1:M)
		v = v * .k (kpdf, bw [ncon + j], x [,ncon + j], u [j])
	.sumk (is.weighted, n, w, v) / constant
}

#used in both c/mvc
.cdfc.cks.eval.scalar = function (constant, v, M, ncon, is.weighted, kcdf, bw, n, x, w, low, constv, u)
{	for (j in 1:M)
		v = v * .K (kcdf, bw [ncon + j], x [,ncon + j], u [j], low [j], constv [,j])
	.sumk (is.weighted, n, w, v) / constant
}

#computes denominator, in each interation
.cdfc4chqf.cks.eval.scalar = function (ncon, is.weighted, conditions, kpdf, kcdf, bw, n, x, w, u)
{	v = .precompute.cksc.v (ncon, conditions, kpdf, bw, n, x)
	btm = .sumk (is.weighted, n, w, v)
	v = v * .K (kcdf, bw [ncon + 1], x [,ncon + 1], u)
	top = .sumk (is.weighted, n, w, v)
	top / btm
}

#calls .qfc.from.chqf.cks -> .cdfc4chqf.cks.eval -> .cdfc4chqf.cks.eval.scalar
.chqf.cks.eval.ext = function (chqf.f, y)
{	n = nrow (y)
	m = chqf.f %$% "m"
	x = matrix (0, n, m)
	x [,1] = (chqf.f %$% "qf1")(y [,1])
	if (m > 1)
	{	for (i in 1:n)
		{	for (j in 2:m)
			{	conditions = x [i, 1:(j - 1)]
				F.inv = .qfc.from.chqf.cks (chqf.f, j, conditions)
				x [i, j] = F.inv (y [i, j])
			}
		}
	}
	x
}

.pwith.eval.uv = function (bw, kcdf, n, x, a, b, isw=FALSE, w)
	.pwith.eval.mv (bw, kcdf, n, 1, cbind (x), a, b, isw, w)

.pwith.eval.mv = function (bw, kcdf, n, m, x, a, b, isw=FALSE, w)
{	y = rep (1, n)
	for (j in 1:m)
	{	area.below.lo = .K (kcdf, bw [j], x [,j], a [j], TRUE)
		area.above.up = .K (kcdf, bw [j], x [,j], b [j], FALSE)
		y = y * (1 - area.below.lo - area.above.up)
	}
	.sumk (isw, n, w, y)
}

.update.wkc.uv = function (ist.lo, bw, kcdf, XLIM, n, x, low)
	as.vector (.update.wkc.mv (ist.lo, bw, kcdf, rbind (XLIM), n, 1, cbind (x), low) )

.update.wkc.mv = function (ist.lo, bw, kcdf, XLIM, n, m, x, low)
{	u = .cdf.lower.side (low, XLIM)
	y = matrix (0, n, m)
	for (j in 1:m)
	{	if (ist.lo [j])
			y [,j] = .K (kcdf, bw [j], x [,j], u [j], low [j])
	}
	y
}

.update.wk.uv = function (ist, bw, kcdf, XLIM, n, x)
	as.vector (.update.wk.mv (rbind (ist), bw, kcdf, rbind (XLIM), n, 1, cbind (x) ) )

.update.wk.mv = function (ist, bw, kcdf, XLIM, n, m, x)
{	a = XLIM [,1]
	b = XLIM [,2]
	w = rep (1, n)
	for (j in 1:m)
	{	area.below.lo = area.above.up = 0
		if (ist [j, 1])
			area.below.lo = .K (kcdf, bw [j], x [,j], a [j], TRUE)
		if (ist [j, 2])
			area.above.up = .K (kcdf, bw [j], x [,j], b [j], FALSE)
		w = w * (1 - area.below.lo - area.above.up)
	}
	w
}
