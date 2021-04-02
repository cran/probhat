#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.binary.seq = function (m, n)
{	rep.vector = round (2 ^ (0:(m - 1) ) )
	verts = matrix (0, nrow=n, ncol=m)
	for (j in 1:m)
		verts [,j] = rep (c (0, 1), times=rep.vector [j], each=rep.vector [1 + m - j])
	verts
}

.binary.sign = function (binary)
{   sums = apply (binary, 1, sum)
    sign = rep (1L, nrow (binary) )
    if (ncol (binary) %% 2 == 0)
        sign [(1 + sums) %% 2 == 0] = -1L
    else
        sign [sums %% 2 == 0] = -1L
    sign
}

probmv = function (sf, a, b)
{	if (inherits (sf, "ccdfmv") )
	{	if (all (sf %$% ".low") )
			.comb.prob (sf, a, b)
		else
			stop ("needs cdf based on lower tails")
	}
	else
		stop ("needs ccdfmv object")
}

.comb.prob = function (F, a, b)
{   has.matrix.args = (is.matrix (a) && is.matrix (b) )
    if (has.matrix.args)
    {   m = ncol (a)
        if (nrow (a) != nrow (b) )
            stop ("nrow (a) needs to equal nrow (b)")
        if (ncol (a) != ncol (b) )
            stop ("ncol (a) needs to equal ncol (b)")
    }
    else
    {   m = length (a)
        if (length (a) != length (b) )
            stop ("length (a) needs to equal length (b)")
    }
    nF = as.integer (round (2 ^ m) )
    binary = .binary.seq (m, nF)
    sign = .binary.sign (binary)
    y = 0
    for (i in 1:nF)
    {   x = a
        j = as.logical (binary [i,])
        if (has.matrix.args)
            x [,j] = b [,j]
        else
            x [j] = b [j]
        y = y + sign [i] * F (x)
    }
    y
}

pwith = function (...)
	UseMethod ("pwith")

pwith.cksuv = function (sf, xlim = c (a, b), ..., a=-Inf, b=Inf)
{	. = attributes (sf)
	xlim = .val.XLIM.uv (xlim)
	data = .$data
	.pwith.eval.uv (.$bw, .$kernel@F, data$n, data$x, xlim [1], xlim [2],
		isw=.$is.weighted, data$w)
}

pwith.cksmv = function (sf, xlim = cbind (a, b), ..., a=-Inf, b=Inf)
{	. = attributes (sf)
	xlim = .val.XLIM.mv (xlim, .$m)
	data = .$data
	.pwith.eval.mv (.$bw, .$kernel@F, data$n, .$m, data$x, xlim [,1], xlim [,2], isw=.$is.weighted, data$w)
}

psv = function (fh)
{	x = attr (fh, "data")$x
	fh (x)
}

ph4.pcomp2 = function (fh, gh, sqrt.mse=TRUE, aggregate=TRUE, dfh = psv (fh), dgh = psv (gh) )
{	nf = length (dfh)
	ng = length (dgh)
	xf = attr (fh, "data")$x
	xg = attr (gh, "data")$x
	kf = sum ( (dfh - gh (xf) )^2) / nf #f data, evaluated with f and g
	kg = sum ( (dgh - fh (xg) )^2) / ng #g data, evaluated with f and g
	if (sqrt.mse)
	{	kf = sqrt (kf)
		kg = sqrt (kg)
	}
	if (aggregate)
		(kf + kg) / 2
	else
		c (kf, kg)
}

pdist = function (fhs, ..., sqrt.mse=TRUE)
{	if (! is.list (fhs) )
		stop ("list required")
	nfh = length (fhs)
	fv0 = vector ("list", nfh)
	dists = matrix (0, nfh, nfh)
	for (i in seq_len (nfh) )
	{	if (! is.pdf (fhs [[i]]) )
			stop ("fhs should only contain (probhat) density functions")
		fv0 [[i]] = psv (fhs [[i]])
	}
	if (nfh > 1)
	{	for (i in 1:(nfh - 1) )
		{	for (j in (i + 1):nfh)
			{	dists [i, j] = dists [j, i] =
					ph4.pcomp2 (fhs [[i]], fhs [[j]], sqrt.mse, TRUE, fv0 [[i]], fv0 [[j]])
			}
		}
	}
	if (inherits (fhs, "ph4.gset") )
		rownames (dists) = colnames (dists) = attr (fhs, "levnames")
	dists
}
