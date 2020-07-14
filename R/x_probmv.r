#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

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
		.comb.prob (sf, a, b)
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
