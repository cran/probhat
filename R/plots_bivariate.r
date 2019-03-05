.contour = function (x, y, z, xlab="x", ylab="y", ...)
	contour (x, y, z, xlab=xlab, ylab=ylab, ...)

plot.nppdfmv = function (x, use.plot3d=FALSE, xlab, ylab, npoints=30, ..., all=FALSE)
{	nppdfmv.f = x

	. = attributes (nppdfmv.f)
	n = .$n
	if (.$m != 2)
		stop ("can only plot nppdfmv if bivariate")
	if (all)
	{	par (mfrow=c (2, 2) )
        p0 = par (mar=c (2, 2.5, 1, 0.175) )
        npcdfmv.f = function (x) {.npcdfmv.eval (x)}
		attributes (npcdfmv.f) = attributes (nppdfmv.f)
        plot (nppdfmv.f, FALSE, xlab="", ylab="", npoints, drawlabels=FALSE, ...)
        plot (nppdfmv.f, TRUE, xlab="", ylab="", npoints, ...)
        plot (npcdfmv.f, FALSE, xlab="", ylab="", npoints, drawlabels=FALSE, ...)
        plot (npcdfmv.f, TRUE, xlab="", ylab="", npoints, ...)
        par (p0)
	}
	else
	{	if (missing (xlab) )
			xlab = .$varnames [1]
		if (missing (ylab) )
			ylab = .$varnames [2]
		xrng = range (.$x [,1]) + c (-0.5, 0.5) * .$bw [1]
		yrng = range (.$x [,2]) + c (-0.5, 0.5) * .$bw [2]
		x = seq (xrng [1], xrng [2], length.out=npoints)
		y = seq (yrng [1], yrng [2], length.out=npoints)
		z = outer (x, y, .mix, nppdfmv.f)
		if (use.plot3d)
			plot3d.surf (x, y, z, xlab=xlab, ylab=ylab, ...)
		else
			.contour (x, y, z, xlab, ylab, ...)
	}
	
}

plot.npcdfmv = function (x, use.plot3d=FALSE, xlab, ylab, npoints=30, ...)
{	npcdfmv.f = x

	. = attributes (npcdfmv.f)
	n = .$n
	if (.$m != 2)
		stop ("can only plot npcdfmv if bivariate")
	if (missing (xlab) )
		xlab = .$varnames [1]
	if (missing (ylab) )
			ylab = .$varnames [2]
	xrng = range (.$x [,1]) + c (-0.5, 0.5) * .$bw [1]
	yrng = range (.$x [,2]) + c (-0.5, 0.5) * .$bw [2]
	x = seq (xrng [1], xrng [2], length.out=npoints)
	y = seq (yrng [1], yrng [2], length.out=npoints)
	z = outer (x, y, .mix, npcdfmv.f)
	if (use.plot3d)
		plot3d.surf (x, y, z, xlab=xlab, ylab=ylab, ...)
	else
		.contour (x, y, z, xlab, ylab, ...)
}

plot.chained.npcdfmv.inverse = function (x, ...)
	stop ("can't plot chained.npcdfmv.inverse")

.mix = function (x, y, f)
	f (cbind (x, y) )
