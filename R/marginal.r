marginal = function (npmv)
{	. = attributes (npmv)
	muv = vector ("list", .$m)
	class (muv) = "marginal"
	attr (muv, "varnames") = .$varnames
	if (inherits (npmv, "nppdfmv") )
	{	for (j in 1:.$m)
			muv [[j]] = nppdfuv (.$x [,j], kernel.pdf=.$kernel.pdf, bw=.$bw [j], w=.$w)
	}
	else if (inherits (npmv, "npcdfmv") )
	{	for (j in 1:.$m)
			muv [[j]] = npcdfuv (.$x [,j], kernel.cdf=.$kernel.cdf, bw=.$bw [j], w=.$w)
	}
	else
		stop ("marginal needs nppdfmv or npcdfmv object")
	muv
}

plot.marginal = function (x, nrow, ncol, with.points=FALSE, ...)
{	muv = x

	m = length (muv)
	varnames = attr (muv, "varnames")
	if (missing (nrow) || missing (ncol) )
		nrow = ncol = ceiling (sqrt (m) )
	p0 = par (mfrow = c (nrow, ncol) )
	for (j in 1:m)
		plot (muv [[j]], with.points, main=varnames [j], ...)
	par (p0)
}
