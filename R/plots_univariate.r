.plot.nested = function (., xrng, x, y, nmax, main="", xlab="x", ylab="y")
{	ymax = max (y)
	ymin = -ymax / 4
	ymax2 = -ymax / 20

	plot.new ()
	plot.window (xlim=xrng, ylim=c (ymin, ymax) )
	title (main=main, xlab=xlab, ylab=ylab)
	box ()
	abline (h=0)
	axis (1)
	axis (2, round (seq (0, ymax, length.out=4), 2) )
	
	lines (x, y)

	n = .$n
	x = .$x
	if (n > nmax)
	{	n = nmax
		x = sample (x, n)
	}
	y = runif (n, ymin, ymax2)
	points (x, y)
}

plot.nppdfuv = function (x, with.points=FALSE, nmax=2000, ...)
{   nppdfuv.f = x

    . = attributes (nppdfuv.f)
	xrng = range (.$x) + c (-0.5, 0.5) * .$bw
    x = seq (xrng [1], xrng [2], length.out=200)
    y = nppdfuv.f (x)
	if (with.points)
		.plot.nested (., xrng, x, y, nmax, ...)
	else
		plot (x, y, type="l", ...)
}

lines.nppdfuv = function (x, ...)
{   nppdfuv.f = x

    . = attributes (nppdfuv.f)
    xrng = range (.$x) + c (-0.5, 0.5) * .$bw
    x = seq (xrng [1], xrng [2], length.out=200)
    y = nppdfuv.f (x)
    lines (x, y, type="l", ...)
}

plot.npcdfuv = function (x, with.points=FALSE, nmax=2000, ...)
{   npcdfuv.f = x

	. = attributes (npcdfuv.f)
    xrng = range (.$x) + c (-0.5, 0.5) * .$bw
    x = seq (xrng [1], xrng [2], length.out=200)
	y = npcdfuv.f (x)
    if (with.points)
		.plot.nested (., xrng, x, y, nmax, ...)
	else
		plot (x, y, type="l", ...)
}

lines.npcdfuv = function (x, ...)
{   npcdfuv.f = x

    . = attributes (npcdfuv.f)
    xrng = range (.$x) + c (-0.5, 0.5) * .$bw
    x = seq (xrng [1], xrng [2], length.out=200)
	y = npcdfuv.f (x)
	lines (x, y, type="l", ...)
}

plot.npcdfuv.inverse = function (x, ...)
{   npcdfuv.f.inverse = x

    . = attributes (npcdfuv.f.inverse)
    y = seq (0, 1, length.out=200)
    x = npcdfuv.f.inverse (y)
    plot (y, x, type="l", ...)
}

lines.npcdfuv.inverse = function (x, ...)
{   npcdfuv.f.inverse = x

    . = attributes (npcdfuv.f.inverse)
    y = seq (0, 1, length.out=200)
    x = npcdfuv.f.inverse (y)
    lines (y, x, type="l", ...)
}
