#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.kpdf.c = function (x, y)
{	y [x < -1 | x > 1] = 0
	y
}

.kcdf.c = function (x, y)
{	y [x < -1] = 0
	y [x > 1] = 1
	y
}

.f2 = function (f1, f2, x)
{	I1 = (x >= -1 & x < 0)
	I2 = (x >= 0 & x <= 1)
	x1 = x [I1]
	x2 = x [I2]
	y = rep (0, length (x) )
	y [I1] = f1 (x1)
	y [I2] = f2 (x2)
	y
}

.dkernel = function (class, name, pmf, cdf, hbw)
{	r = c (-hbw, hbw)
	pmf = EXTEND (pmf, .KV.pmf, hbw, xlim=r, freq=FALSE)
	cdf = EXTEND (cdf, .KV.dcdf, hbw, xlim=r)
	EXTEND (LIST (pmf, cdf), c (class, "dkernel", "kernel"), name)
}

.ckernel = function (class, name, pdf, cdf)
{	r = c (-1, 1)
	pdf = EXTEND (pdf, .KV.pdf, xlim=r)
	cdf = EXTEND (cdf, .KV.ccdf, xlim=r)
	EXTEND (LIST (pdf, cdf), c (class, "ckernel", "kernel"), name)
}

binomial.dkernel = function (bw)
{	n = bw - 1
	if (n %% 2 == 0)
	{	hbw = n %/% 2
		.dkernel ("binomial.dkernel", "binomial", .binomial.pmf, .binomial.cdf, hbw)
	}
	else stop ("unsuitable bandwidth")
}

uniform.ckernel = function ()
	.ckernel ("uniform.ckernel", "uniform", .uniform.pdf, .uniform.cdf)
triangular.ckernel = function ()
	.ckernel ("triangular.ckernel", "triangular", .triangular.pdf, .triangular.cdf)
epanechnikov.ckernel = function ()
	.ckernel ("epanechnikov.ckernel", "epanechnikov", .epanechnikov.pdf, .epanechnikov.cdf)
truncnorm.ckernel = function ()
	.ckernel ("truncnormal.ckernel", "truncated normal", .truncnorm.pdf, .truncnorm.cdf)
biweight.ckernel = function ()
	.ckernel ("biweight.kernel", "biweight", .biweight.pdf, .biweight.cdf)
triweight.ckernel = function ()
	.ckernel ("triweight.ckernel", "triweight", .triweight.pdf, .triweight.cdf)
tricube.ckernel = function ()
	.ckernel ("tricube.ckernel", "tricube", .tricube.pdf, .tricube.cdf)
bell.spline = function ()
	.ckernel ("bell.spline", "bell spline", .bell.spline.pdf, .bell.spline.cdf)

.binomial.pmf = function (x, ...)
{	. = THAT ()
	with (., dbinom (x + hbw, 2 * hbw, 0.5) )
}

.binomial.cdf = function (x)
{	. = THAT ()
	with (., pbinom (x + hbw, 2 * hbw, 0.5) )
}

.uniform.pdf = function (x)
{	y = rep (0.5, length (x) )
	.kpdf.c (x, y)
}

.uniform.cdf = function (x)
{	y = 0.5 + 0.5 * x
	.kpdf.c (x, y)
}

.triangular.pdf = function (x)
{	y = 1 - abs (x)
	.kpdf.c (x, y)
}

.triangular.cdf = function (x)
	.f2 (.triangular.cdf.1, .triangular.cdf.2, x)
.triangular.cdf.1 = function (x)
	0.5 + x + 0.5 * x ^ 2
.triangular.cdf.2 = function (x)
	0.5 + x - 0.5 * x ^ 2

.epanechnikov.pdf = function (x)
{	y = 0.75 * (1 - x ^ 2)
	.kpdf.c (x, y)
}

.epanechnikov.cdf = function (x)
{	y = 0.75 * (x - 1 / 3 * x ^ 3) + 0.5
	.kcdf.c (x, y)
}

.truncnorm.pdf = function (x)
{	y = 2.8070337683438038 * ( (dnorm (2.8070337683438038 * x) - 0.0077608934080090) / 0.995)
	.kpdf.c (x, y)
}
.truncnorm.cdf = function (x)
{	y = (pnorm (2.8070337683438038 * x) - 0.0025) / 0.995
	.kcdf.c (x, y)
}

.biweight.pdf = function (x)
{	y = 0.9375 * (1 - x ^ 2) ^ 2
	.kpdf.c (x, y)
}

.biweight.cdf = function (x)
{	y = 0.9375 * (x - 2 / 3 * x ^ 3 + 0.2 * x ^ 5) + 0.5
	.kcdf.c (x, y)
}

.triweight.pdf = function (x)
{	 y = 1.09375 * (1 - x ^ 2) ^ 3
	.kpdf.c (x, y)
}

.triweight.cdf = function (x)
{	y = 1.09375 * (x - x ^ 3 + 0.6 * x ^ 5 - 1 / 7 * x ^ 7) + 0.5
	.kcdf.c (x, y)
}

.tricube.pdf = function (x)
{	y = 0.8641975 * (1 - abs (x) ^ 3) ^ 3
	.kpdf.c (x, y)
}

.tricube.cdf = function (x)
	.f2 (.tricube.cdf.1, .tricube.cdf.2, x)
.tricube.cdf.1 = function (x)
	0 + 70 / 81 * (x + 0.75 * x ^ 4 + 3 / 7 * x ^ 7 + 0.1 * x ^ 10) + 0.5
.tricube.cdf.2 = function (x)
	0 + 70 / 81 * (x - 0.75 * x ^ 4 + 3 / 7 * x ^ 7 - 0.1 * x ^ 10) + 0.5

.bell.spline.pdf = function (x)
{	y = rep (0, length (x) )
	#subinterval 1
	I = (x > -1 & x < -0.5)
	y [I] = 2 + 4 * x [I] + 2 * x [I] ^ 2
	#center and subintervals 2 and 3
	I = (x >= -0.5 & x <= 0.5)
	y [I] = 1 - 2 * x [I] ^ 2
	#subinterval 4
	I = (x > 0.5 & x < 1)
	y [I] = 2 - 4 * x [I] + 2 * x [I] ^ 2
	y
}

.bell.spline.cdf = function (x)
{	y = rep (0, length (x) )
	#subinterval 1
	I = (x > -1 & x < -0.5)
	y [I] = 2 / 3 + 2 * x [I] + 2 * x [I] ^ 2 + 2 / 3 * x [I] ^ 3
	#center and subintervals 2 and 3
	I = (x >= -0.5 & x <= 0.5)
	y [I] = 0.5 + x [I] - 2 / 3 * x [I] ^ 3
	#subinterval 4
	I = (x > 0.5 & x < 1)
	y [I] = 1 / 3 + 2 * x [I] - 2 * x [I] ^ 2 + 2 / 3 * x [I] ^ 3
	#x >= 1
	I = (x >= 1)
	y [I] = 1
	y
}
