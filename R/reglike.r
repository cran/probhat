.reglike.varindex = function (varname, varnames)
	match (toupper (varname), toupper (varnames) )

.reglike = function (model, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, nc, nc.npc, smoothness, bw, w, x)
{	.test.nc (nc)
	model = as.list (model)
	stat.exp = as.character (model [[2]])
	stat = stat.exp [1]
	reglike.y.varname = stat.exp [2]
	reglike.x.varname = as.character (model [[3]])

	m = ncol (x)
	varnames = colnames (x)
	if (is.null (varnames) || length (unique (toupper (varnames) ) ) != m) 
		stop ("x needs unique column names")
	Jx = .reglike.varindex (reglike.x.varname, varnames)
	Jy = .reglike.varindex (reglike.y.varname, varnames)
	if (is.na (Jx) || is.na (Jy) )
		stop ("variables in formula not in data")
	J = 1:m
	J = c (Jx, J [-c (Jx, Jy)], Jy)

	if (!missing (smoothness) && length (smoothness) == m)
		smoothness = smoothness [J]
	if (!missing (bw) && length (bw) == m)
		bw = bw [J]
	. = .npmv (smoothness, bw, w, x [,J])

	conditions = numeric (m - 1)
	if (m - 1 > 1)
	{	for (j in 2:(m - 1) )
			conditions [j] = eccdf.inverse (.$x [,j], w)(0.5)
	}

	xrng = range (.$x [,1])
	cx = seq (xrng [1], xrng [2], length.out=nc)
	cy = numeric (nc)

	if (stat [1] == "mean")
	{	for (i in 1:nc)
		{	conditions [1] = cx [i]
			cy [i] = npmean (nppdfc (conditions, .$x,, kernel.pdf, kernel.cdf, nc.npc, bw=.$bw, w=w) )
		}
	}
	else if (stat [1] == "mode")
	{	for (i in 1:nc)
		{	conditions [1] = cx [i]
			cy [i] = npmode (nppdfc (conditions, .$x,, kernel.pdf, kernel.cdf, nc.npc, bw=.$bw, w=w) )
		}
	}
	else if (stat [1] == "median")
	{	for (i in 1:nc)
		{	conditions [1] = cx [i]
			cy [i] = npcdfc.inverse (conditions, .$x, kernel.pdf, kernel.cdf, nc.npc, bw=.$bw, w=w)(0.5)
		}
	}
	else if (stat [1] == "quantile")
	{	for (i in 1:nc)
		{	p = as.numeric (stat.exp [3])
			conditions [1] = cx [i]
			cy [i] = npcdfc.inverse (conditions, .$x, kernel.pdf, kernel.cdf, nc.npc, bw=.$bw, w=w)(p)
		}
	}
	else
		stop ("unsuitable stat argument")
	conditions [1] = NA

	.$conditions = conditions
	.$cx = cx
	.$cy = cy
	if (any (is.na (cy) ) )
		stop ("\nconditional distribution evaluates to NA\n(possibly denominator with zero density)")
	.$ct = .slopes (nc, cx, cy)

	.
}

reglike = function (model, x, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, nc=30, nc.npc=nc, smoothness=0.65, bw, w=NA)
{	.= .reglike (model, kernel.pdf, kernel.cdf, nc, nc.npc, smoothness, bw, w, x)
	reglike.f = function (x) {.reglike.eval (x)}
	attributes (reglike.f) = list (
		class = "reglike",
		model=model, conditions=.$conditions,
		spline=TRUE, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=nc, nc.npc=nc.npc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
	reglike.f
}

.reglike.eval = function (x)
{	. = attributes (sys.function (-1) )
	n = length (x)
    y = numeric (n)
    for (i in 1:n)
		y [i] = .spline.eval (.$nc, .$cx, .$cy, .$ct, x [i])
    y
}

plot.reglike = function (x, with.points=FALSE, ...)
{	reglike.f = x

	. = attributes (reglike.f)
	if (with.points)
	{	x = .$x [,1]
		y = .$x [,.$m]
		plot (x, y, ...)
		lines (.$cx, .$cy)
	}
	else
	{	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
		y = reglike.f (x)
		plot (x, y, type="l", ...)
	}
}

lines.reglike = function (x, ...)
{	reglike.f = x

	. = attributes (reglike.f)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	y = reglike.f (x)
	lines (x, y, ...)
}
