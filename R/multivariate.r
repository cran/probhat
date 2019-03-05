.npmv = function (smoothness, bw, w, x)
{   if (!is.matrix (x) )
		stop ("multivariate models need matrix")
	n = nrow (x)
    m = ncol (x)
	varnames = colnames (x)
	if (is.null (varnames) || length (unique (toupper (varnames) ) ) != m) 
		stop ("x needs unique column names")
	af = all (is.finite (x) )
	if (!af)
        stop ("no missing values allowed")
	for (j in 1:m)
		if (length (unique (x [,j]) ) < 2)
			stop ("each column of x needs\n2 or more unique values")

	weighted = (!is.na (w [1]) )
	if (weighted)
	{	w = as.numeric (w)
		if (n != length (w) )
			stop ("length(w) not equal nrow(x)")
		if (any (is.na (w) ) )
			stop ("no missing values allowed")
		if (round (sum (w), 1) != 1)
			stop ("sum(w) must be approx 1")
		if (any (w <= 0) )
			stop ("all w need to be > 0")
	}

	if (!missing (smoothness) )
	{	m.2 = length (smoothness)
		if (m.2 == 1)
			smoothness = rep (smoothness, m)
		else if (m.2 != m)
			stop ("length(smoothness) must equal 1 or m")
	}
	if (missing (bw) )
	{	bw = numeric (m)
		for (j in 1:m)
			bw [j] = smoothness [j] * diff (range (x [,j]) )
	}
	else
	{	smoothness = NA
		if (m != length (bw) )
			stop ("length(bw) must equal m")
	}

	list (weighted=weighted, smoothness=smoothness, bw=bw, n=n, m=m, varnames=varnames, w=w, x=x)
}

nppdfmv = function (x, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, smoothness=0.65, bw, w=NA)
{	. = .npmv (smoothness, bw, w, x)
	nppdfmv.f = function (x) {.nppdfmv.eval (x)}
	attributes (nppdfmv.f) = list (
		class="nppdfmv",
		kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
	if (.$m == 1)
		stop ("multivariate models need 2 or more variables")
    nppdfmv.f
}

.nppdfmv.eval = function (x)
{   . = attributes (sys.function (-1) )
    if (!is.matrix (x) )
		x = rbind (x)
    n = nrow (x)
    m = ncol (x)
    if (.$m != m)
        stop ("incorrect number of columns")
	y = numeric (n)
	for (i in 1:n)
		y [i] = .nppdfmv.eval.ext (., x [i,])
	y
}

npcdfmv = function (x, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, smoothness=0.65, bw, w=NA)
{   . = .npmv (smoothness, bw, w, x)
	npcdfmv.f = function (x) {.npcdfmv.eval (x)}
    attributes (npcdfmv.f) = list (
		class="npcdfmv",
		kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
	if (.$m == 1)
		stop ("multivariate models need 2 or more variables")
    npcdfmv.f
}

.npcdfmv.eval = function (x)
{	. = attributes (sys.function (-1) )
	if (!is.matrix (x) )
		x = rbind (x)
    n = nrow (x)
	m = ncol (x)
    if (.$m != m)
        stop ("incorrect number of columns")
    y = numeric (n)
	for (i in 1:n)
		y [i] = .npcdfmv.eval.ext (., x [i,])
	y
}

chained.npcdfmv.inverse = function (x, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, smoothness=0.65, bw, w=NA)
{   . = .npmv (smoothness, bw, w, x)
	chained.f = function (x) {.chained.npcdfmv.inverse.eval (x)}
    attributes (chained.f) = list (
		class="chained.npcdfmv.inverse",
		kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
	if (.$m == 1)
		stop ("multivariate models need 2 or more variables")
    chained.f
}

.chained.npcdfmv.inverse.eval = function (y)
{	. = attributes (sys.function (-1) )
	if (!is.matrix (y) )
		y = rbind (y)
	n = nrow (y)
	m = ncol (y)
    if (.$m != m)
        stop ("incorrect number of columns")
	x = matrix (0, nrow=n, ncol=m)
	colnames (x) = .$varnames
	x [,1] = npcdfuv.inverse (.$x [,1], kernel.pdf=.$kernel.pdf, bw=.$bw [1], w=.$w)(y [,1])
	for (i in 1:n)
	{	for (j in 2:m)
		{	conditions = x [i, 1:(j - 1)]
			x [i, j] = npcdfc.inverse (conditions, .$x [,1:j], kernel.pdf=.$kernel.pdf, bw=.$bw [1:j], w=.$w)(y [i, j])
		}
	}
	x
}

.nppdfmv.eval.ext = function (., x)
{	y = 0
	for (i in 1:.$n)
	{	y.local = 1
		for (j in 1:.$m)
		{	dist = x [j] - .$x [i, j]
			l = 2 / .$bw [j] * .$kernel.pdf (2 / .$bw [j] * dist)
			y.local = y.local * l
		}
		if (.$weighted)
			y = y + .$w [i] * y.local
		else
			y = y + y.local
	}
	if (.$weighted)
		y
	else
		y / .$n
}

.npcdfmv.eval.ext = function (., x)
{	y = 0
	for (i in 1:.$n)
	{	y.local = 1
		for (j in 1:.$m)
		{	dist = x [j] - .$x [i, j]
			L = .$kernel.cdf (2 / .$bw [j] * dist)
			y.local = y.local * L
		}
		if (.$weighted)
			y = y + .$w [i] * y.local
		else
			y = y + y.local
	}
	if (.$weighted)
		y
	else
		y / .$n
}

.emixed.eval.ext = function (., x)
{	y = 0
	for (i in 1:.$n)
	{	y.local = 1
		for (j in 1:(.$m - 1) )
		{	dist = x [j] - .$x [i, j]
			l = 2 / .$bw [j] * .$kernel.pdf (2 / .$bw [j] * dist)
			y.local = y.local * l
		}
		dist = x [.$m] - .$x [i, .$m]
		L = .$kernel.cdf (2 / .$bw [.$m] * dist)
		y.local = y.local * L

		if (.$weighted)
			y = y + .$w [i] * y.local
		else
			y = y + y.local
	}
	if (.$weighted)
		y
	else
		y / .$n
}
