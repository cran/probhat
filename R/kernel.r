#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.KV.pmf = c ("pmfuv", "dpduv", "pmf", "dpd", "phpd")
.KV.dcdf = c ("dcdfuv", "dpduv", "dcdf", "dpd", "phpd")
.KV.pdf = c ("pdfuv", "cpduv", "pdf", "cpd", "phpd")
.KV.ccdf = c ("ccdfuv", "cpduv", "ccdf", "cpd", "phpd")

setOldClass ("phob")
setOldClass (c (.KV.pmf, "function", "phob") )
setOldClass (c (.KV.dcdf, "function", "phob") )
setOldClass (c (.KV.pdf, "function", "phob") )
setOldClass (c (.KV.ccdf, "function", "phob") )

setClass ("Kernel", contains = c ("VIRTUAL", "phob"),
	slots = list (name="character", f="function", F="function") )
setClass ("DKernel", contains = c ("Kernel"),
	slots = list (xlim="numeric") )
setClass ("CKernel", contains = c ("Kernel") )

setClass ("Uniform.CKernel", contains="CKernel")
setClass ("Triangular.CKernel", contains="CKernel")
setClass ("Epanechnikov.CKernel", contains="CKernel")
setClass ("TrGaussian.CKernel", contains="CKernel")
setClass ("Biweight.CKernel", contains="CKernel")
setClass ("Triweight.CKernel", contains="CKernel")
setClass ("Tricube.CKernel", contains="CKernel")
setClass ("Bell.Spline", contains="CKernel")

setMethod ("show", "Kernel", function (object) print (object) )

discretized.kernel = function (n, ck=BIWEIGHT.CKERNEL, ..., xlim)
{	name = paste ("Discretized", ck@name)

	f = function (x, ...)
	{	nx = length (x)
		if (nx == 0)
			y = numeric ()
		else
		{	I = (x >= xlim [1] & x <= xlim [2])
			y = rep (0, nx)
			y [I] = .probs [x [I] - xlim [1] + 1]
		}
		y
	}
	
	F = function (x, ...)
	{	nx = length (x)
		if (nx == 0)
			y = numeric ()
		else
		{	I = (x >= xlim [1] & x <= xlim [2])
			U = (x > xlim [2])
			y = rep (0, nx)
			y [I] = .PROBS [x - xlim [1] + 1]
			y [U] = 1
		}
		y
	}

	f = .EXTEND (f, .KV.pmf)
	F = .EXTEND (F, .KV.dcdf)

	if (missing (xlim) )
	{	n = as.integer (n)
		if (length (n) != 1)
			stop ("length (n) != 1")
		if (n <= 0 || n %% 2 == 0)
			stop ("n needs to be a odd positive integer, or provide xlim")
		hn = (n - 1L) %/% 2
		xlim = c (-hn, hn)
	}
	else
	{	xlim = as.integer (xlim)
		if (length (xlim) != 2)
			stop ("length (xlim) != 2")
		if (diff (xlim) <= 0)
			stop ("currently, needs increasing xlim values")
		n = xlim [2] - xlim [1] + 1
	}

	u = .midpoints (seq (-1, 1, length.out = n + 1) )
	.probs = ck@f (u)
	.probs = .probs / sum (.probs)
	.PROBS = cumsum (.probs)

	new ("DKernel", name=name, f=f, F=F, xlim=xlim)
}

.con.ck = function (name, pdf, cdf, ..., con = paste0 (name, ".CKernel") ) 
{	f = function (x)
	{	y = pdf (x)
		y [x < -1 | x > 1] = 0
		y
	}
	F = function (x)
	{	y = cdf (x)
		y [x < -1] = 0
		y [x > 1] = 1
		y
	}

	f = .EXTEND (f, .KV.pdf)
	F = .EXTEND (F, .KV.ccdf)

	new (con, name=name, f=f, F=F)
}

.F2 = function (F1, F2, x)
{	I1 = (x >= -1 & x < 0)
	I2 = (x >= 0 & x <= 1)
	x1 = x [I1]
	x2 = x [I2]
	y = rep (0, length (x) )
	y [I1] = F1 (x1)
	y [I2] = F2 (x2)
	y
}

ph.printf.Kernel = function (k, ...)
	cat (k@name, "Kernel\n")

.uniform.pdf = function (x) rep (0.5, length (x) )
.uniform.cdf = function (x) 0.5 + 0.5 * x

.triangular.pdf = function (x) 1 - abs (x)
.triangular.cdf = function (x) .F2 (function (x) 0.5 + x + 0.5 * x ^ 2, function (x) 0.5 + x - 0.5 * x ^ 2, x)

.epanechnikov.pdf = function (x) 0.75 * (1 - x ^ 2)
.epanechnikov.cdf = function (x) 0.75 * (x - 1 / 3 * x ^ 3) + 0.5

.truncnorm.pdf = function (x) 2.8070337683438038 * ( (dnorm (2.8070337683438038 * x) - 0.0077608934080090) / 0.995)
.truncnorm.cdf = function (x) (pnorm (2.8070337683438038 * x) - 0.0025) / 0.995

.biweight.pdf = function (x) 0.9375 * (1 - x ^ 2) ^ 2
.biweight.cdf = function (x) 0.9375 * (x - 2 / 3 * x ^ 3 + 0.2 * x ^ 5) + 0.5

.triweight.pdf = function (x) 1.09375 * (1 - x ^ 2) ^ 3
.triweight.cdf = function (x) 1.09375 * (x - x ^ 3 + 0.6 * x ^ 5 - 1 / 7 * x ^ 7) + 0.5

.tricube.pdf = function (x) 0.8641975 * (1 - abs (x) ^ 3) ^ 3
.tricube.cdf = function (x)
{	.F2 (
		function (x) 0 + 70 / 81 * (x + 0.75 * x ^ 4 + 3 / 7 * x ^ 7 + 0.1 * x ^ 10) + 0.5,
		function (x) 0 + 70 / 81 * (x - 0.75 * x ^ 4 + 3 / 7 * x ^ 7 - 0.1 * x ^ 10) + 0.5,
		x)
}

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

UNIFORM.CKERNEL = .con.ck ("Uniform", .uniform.pdf, .uniform.cdf)
TRIANGULAR.CKERNEL = .con.ck ("Triangular", .triangular.pdf, .triangular.cdf)
EPANECHNIKOV.CKERNEL = .con.ck ("Epanechnikov", .epanechnikov.pdf, .epanechnikov.cdf)
TRGAUSSIAN.CKERNEL = .con.ck ("Truncated Gaussian", .truncnorm.pdf, .truncnorm.cdf,
	con="TrGaussian.CKernel")
BIWEIGHT.CKERNEL = .con.ck ("Biweight", .biweight.pdf, .biweight.cdf)
TRIWEIGHT.CKERNEL = .con.ck ("Triweight", .triweight.pdf, .triweight.cdf)
TRICUBE.CKERNEL = .con.ck ("Tricube", .tricube.pdf, .tricube.cdf)

BELL.SPLINE = new ("Bell.Spline", name="Bell Spline",
	f = .EXTEND (.bell.spline.pdf, .KV.pdf),
	F = .EXTEND (.bell.spline.cdf, .KV.ccdf) )
