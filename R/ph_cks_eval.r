#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

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

.K = function (kcdf, bw, x, u)
{	dist = u - x
	kcdf (2 / bw * dist)
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

.cdfuv.cks.eval.scalar = function (is.weighted, kcdf, bw, n, x, w, u)
{	y = .K (kcdf, bw, x, u)
	.sumk (is.weighted, n, w, y)
}

.pdfmv.cks.eval.scalar = function (is.weighted, kpdf, bw, n, m, x, w, u)
{	y = rep (1, n)
	for (j in 1:m)
		y = y * .k (kpdf, bw [j], x [,j], u [j])
	.sumk (is.weighted, n, w, y)
}

.cdfmv.cks.eval.scalar = function (is.weighted, kcdf, bw, n, m, x, w, u)
{	y = rep (1, n)
	for (j in 1:m)
		y = y * .K (kcdf, bw [j], x [,j], u [j])
	.sumk (is.weighted, n, w, y)
}

.precompute.cksc.v = function (ncon, conditions, kpdf, bw, n, x)
{	v = rep (1, n)
	for (j in seq_len (ncon) )
		v = v * .k (kpdf, bw [j], x [,j], conditions [j])
	v
}

.pdfc.cks.eval.scalar = function (constant, v, M, ncon, is.weighted, kpdf, bw, n, x, w, u)
{	for (j in 1:M)
		v = v * .k (kpdf, bw [ncon + j], x [,ncon + j], u [j])
	.sumk (is.weighted, n, w, v) / constant
}

.cdfc.cks.eval.scalar = function (constant, v, M, ncon, is.weighted, kcdf, bw, n, x, w, u)
{	for (j in 1:M)
		v = v * .K (kcdf, bw [ncon + j], x [,ncon + j], u [j])
	.sumk (is.weighted, n, w, v) / constant
}

.cdfc.cks.eval.scalar.2 = function (ncon, is.weighted, conditions, kpdf, kcdf, bw, n, x, w, u)
{	v = .precompute.cksc.v (ncon, conditions, kpdf, bw, n, x)
	btm = .sumk (is.weighted, n, w, v)
	v = v * .K (kcdf, bw [ncon + 1], x [,ncon + 1], u)
	top = .sumk (is.weighted, n, w, v)
	top / btm
}

.chqf.cks.eval = function (chqf.f, y)
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
