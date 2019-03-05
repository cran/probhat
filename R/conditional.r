.npc.varindices = function (rv, condnames, varnames)
{	if (length (rv) != 1)
		stop ("length (rv) must equal 1")
	if (length (condnames) != length (varnames) - 1)
		stop ("length (condnames) must be length (varnames) - 1")

	rv = toupper (rv)
	condnames= toupper (condnames)
	varnames = toupper (varnames)

	if (any (rv == condnames) )
		stop ("rv can not be the same as conditioning variables")
	if (length (unique (condnames) ) != length (condnames) )
		stop ("all names must be unique")

	J1 = match (rv, varnames)
	J2 = match (condnames, varnames)
	if (is.na (J1) || any (is.na (J2) ) )
		stop ("variables in formula not in data")

	c (J2, J1)
}

.npc = function (spline, ext, model, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
{	.test.nc (nc)
	if (is.vector (model) )
	{	. = .npmv (smoothness, bw, w, x)
		conditions = as.numeric (model)
		names (conditions) = .$varnames [1:(.$m - 1)]
	}
	else if (inherits (model, "formula") )
	{	model = as.list (model)
		rv = as.character (model [[2]])
		conditions = eval (model [[3]])
		condnames = names (conditions)
		if (is.null (condnames) || length (unique (condnames) ) != length (condnames) ) 
			stop ("conditions need unique names")
		varnames = colnames (x)
		if (is.null (varnames) || length (unique (varnames) ) != length (varnames) ) 
			stop ("x needs unique column names")
		J = .npc.varindices (rv, condnames, varnames)

		m = ncol (x)
		if (!missing (smoothness) && length (smoothness) == m)
			smoothness = smoothness [J]
		if (!missing (bw) && length (bw) == m)
			bw = bw [J]
		. = .npmv (smoothness, bw, w, x [,J])
	}
	else
		stop ("model needs to be vector or formula")

	.$conditions = conditions
	.$nc = .$cx = .$cy = .$ct = NA
	if (.$m < 2)
		stop ("conditional models need m > 1")
	if (spline)
	{	.$conditions = conditions
		.$kernel.pdf = kernel.pdf
		.$kernel.cdf = kernel.cdf

		.$nc = nc
		xrng = range (.$x [,.$m]) + c (-0.5, 0.5) * .$bw [.$m]
		.$cx = seq (xrng [1], xrng [2], length.out=nc)
		.$cy = numeric (nc)
		y = ext (., .$cx [1])
		if (is.na (y) )
			.$cy = .$ct = rep (NA, nc)
		else
		{	.$cy [1] = y
			for (i in 2:nc)
				.$cy [i] = ext (., .$cx [i])
			.$ct = .slopes (nc, .$cx, .$cy)
		}
	}

	.
}

nppdfc = function (model, x, spline=TRUE, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, nc=30, smoothness=0.65, bw, w=NA)
{  	. = .npc (spline, .nppdfc.eval.ext, model, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
	nppdfc.f = function (x) {.nppdfc.eval (x)}
	attributes (nppdfc.f) = list (
		class=c ("nppdfc", "nppdfuv"),
		model=model, conditions=.$conditions,
		spline=spline, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
    nppdfc.f
}

.nppdfc.eval = function (x)
{	. = attributes (sys.function (-1) )
	if (.$spline)
		.npuv.eval.spline (., x, c (0, 0) )
	else
		.npuv.eval.direct (.nppdfc.eval.ext, ., x)
}

npcdfc = function (model, x, spline=TRUE, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, nc=30, smoothness=0.65, bw, w=NA)
{	. = .npc (spline, .npcdfc.eval.ext, model, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
	npcdfc.f = function (x) {.npcdfc.eval (x)}
		attributes (npcdfc.f) = list (
		class=c ("npcdfc", "npcdfuv"),
		model=model, conditions=.$conditions,
		spline=spline, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
    npcdfc.f
}

.npcdfc.eval = function (x)
{	. = attributes (sys.function (-1) )
	if (.$spline)
		.npuv.eval.spline (., x, c (0, 1) )
	else
		.npuv.eval.direct (.npcdfc.eval.ext, ., x)
}

npcdfc.inverse = function (model, x, kernel.pdf=sbc.pdf, kernel.cdf=sbc.cdf, nc=30, smoothness=0.65, bw, w=NA)
{   . = .npc (TRUE, .npcdfc.eval.ext, model, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
	npcdfc.f.inverse = function (y) {.npcdfc.inverse.eval (y)}
	attributes (npcdfc.f.inverse) = list (
		class=c ("npcdfc.inverse", "npcdfuv.inverse"),
		model=model, conditions=.$conditions,
		spline=TRUE, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, m=.$m, varnames=.$varnames, w=.$w, x=.$x)
    npcdfc.f.inverse
}

.npcdfc.inverse.eval = function (y)
{	. = attributes (sys.function (-1) )
	.npcdfuv.inverse.eval.spline (., y)
}

.nppdfc.eval.ext = function (., x)
{	top = .nppdfmv.eval.ext (., c (.$conditions, x) )
	.$m = .$m - 1
	btm = .nppdfmv.eval.ext (., .$conditions)
	if (btm == 0)
		NA
	else
		top / btm
}

.npcdfc.eval.ext = function (., x)
{	top = .emixed.eval.ext (., c (.$conditions, x) )
	.$m = .$m - 1
	btm = .nppdfmv.eval.ext (., .$conditions)
	if (btm == 0)
		NA
	else
		top / btm
}
