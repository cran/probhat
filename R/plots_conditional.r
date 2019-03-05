.xc = function (npc, inverse=FALSE)
{	. = attributes (npc)
	if (.$spline)
	{	if (inverse)
			xrng = range (.$cy)
		else
			xrng = range (.$cx)
	}
	else
		xrng = range (.$x [,.$m]) + c (-0.5, 0.5) * .$bw [.$m]
	seq (xrng [1], xrng [2], length.out=200)
}

plot.nppdfc = function (x, ...)
{	nppdfc.f = x

	x =.xc (nppdfc.f)
    y = nppdfc.f (x)
	plot (x, y, type="l", ...)
}

lines.nppdfc = function (x, ...)
{	nppdfc.f = x

	x =.xc (nppdfc.f)
    y = nppdfc.f (x)
	lines (x, y, ...)
}

plot.npcdfc = function (x, ...)
{	npcdfc.f = x

	x =.xc (npcdfc.f)
    y = npcdfc.f (x)
	plot (x, y, type="l", ...)
}

lines.npcdfc = function (x, ...)
{	npcdfc.f = x

	x =.xc (npcdfc.f)
    y = npcdfc.f (x)
	lines (x, y, ...)
}

plot.npcdfc.inverse = function (x, ...)
{	npcdfc.f.inverse = x

	y = .xc (npcdfc.f.inverse, TRUE)
    x = npcdfc.f.inverse (y)
	plot (y, x, type="l", ...)
}

lines.npcdfc.inverse = function (x, ...)
{	npcdfc.f.inverse = x

	y =.xc (npcdfc.f.inverse, TRUE)
    x = npcdfc.f.inverse (y)
	lines (y, x, ...)
}

