#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.pmfuv.dks.eval = function (x, ..., freq=FALSE, n)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc, .$XLIM)
	y = .iterate.uv (.pmfuv.dks.eval.ext, .$xlim, .$.probs, u=x)
	.scale.freq (y, freq, .$.xsum, n)
}

.cdfuv.dks.eval = function (x, ..., freq=FALSE, n)
{	. = .THAT ()
	x = .val.u.uv (x, .$.any.trunc, .$.is.trunc, .$XLIM)
	y = .iterate.uv (.cdfuv.dks.eval.ext, .$xlim, .$.PROBS, u=x)
	.scale.freq (y, freq, .$.xsum, n)
}

.qfuv.dks.eval = function (p)
{	. = .THAT ()
	.test.y.ok (p)
	.iterate.uv (.qfuv.dks.eval.ext, .$xlim, .$.PROBS, .$.low, u=p)
}

.pmfuv.dks.eval.0 = function (kpmf, x, w, u)
	sum (w * kpmf (u - x) )

.pmfuv.dks.eval.ext = function (xlim, probs, x)
{	if (x < xlim [1] || x > xlim [2])
		0
	else
		probs [x - xlim [1] + 1]
}

.cdfuv.dks.eval.ext = function (xlim, PROBS, x)
{	if (x < xlim [1])
		0
	else if (x > xlim [2])
		1
	else
		PROBS [x - xlim [1] + 1]
}

.qfuv.dks.eval.ext = function (xlim, PROBS, lower=TRUE, y)
{	if (lower)
		q = which (y <= PROBS)[1]
	else
	{	n = length (PROBS)
		q = which (y >= PROBS)[n]
	}
	xlim [1] + q - 1
}
