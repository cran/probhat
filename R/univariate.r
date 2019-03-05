.npuv = function (spline, ext, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x, inverse=FALSE)
{   .test.nc (nc)
	x = as.numeric (x)
	n = length (x)
	af = all (is.finite (x) )
	if (!af)
        stop ("no missing values allowed")
	if (length (unique (x) ) < 2)
		stop ("x needs 2 or more unique values")

	weighted = !is.na (w [1])
	if (weighted)
	{	w = as.numeric (w)
		if (n != length (w) )
			stop ("length(x) not equal length(w)")
		if (any (is.na (w) ) )
			stop ("no missing values allowed")
		if (round (sum (w), 2) != 1)
			stop ("sum(w) must be approx 1")
		if (any (w <= 0) )
			stop ("all w need to be > 0")
	}
	
	rng = range (x)
	if (missing (bw) )
		bw = smoothness * diff (rng)
	else
		smoothness = NA

	. = list (kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=weighted, nc=nc, cx=NA, cy=NA, ct=NA, smoothness=smoothness, bw=bw, n=n, w=w, x=x)
	
	if (spline)
	{	.$cx = seq (rng [1] - bw / 2, rng [2] + bw / 2, length.out=nc)
		.$cy = .npuv.eval.direct (ext, ., .$cx)
		if (inverse)
			.$ct = .slopes (nc, .$cy, .$cx)
		else
			.$ct = .slopes (nc, .$cx, .$cy)
	}
	else
		.$nc = NA

	.
}

nppdfuv = function (x, spline=TRUE, kernel.pdf=sbc.pdf, kernel.cdf=NA, nc=30, smoothness=0.65, bw, w=NA)
{   . = .npuv (spline, .nppdfuv.eval.ext, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
	nppdfuv.f = function (x) {.nppdfuv.eval (x)}
	attributes (nppdfuv.f) = list (
		class="nppdfuv",
		spline=spline, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, w=.$w, x=.$x)
    nppdfuv.f
}

.nppdfuv.eval = function (x)
{	. = attributes (sys.function (-1) )
	if (.$spline)
		.npuv.eval.spline (., x, c (0, 0) )
	else
		.npuv.eval.direct (.nppdfuv.eval.ext, ., x)
}

npcdfuv = function (x, spline=TRUE, kernel.pdf=NA, kernel.cdf=sbc.cdf, nc=30, smoothness=0.65, bw, w=NA)
{   . = .npuv (spline, .npcdfuv.eval.ext, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x)
	npcdfuv.f = function (x) {.npcdfuv.eval (x)}
	attributes (npcdfuv.f) = list (
		class="npcdfuv",
		spline=spline, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cx=.$cx, cy=.$cy, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, w=.$w, x=.$x)
	npcdfuv.f
}

.npcdfuv.eval = function (x)
{	. = attributes (sys.function (-1) )
	if (.$spline)
		.npuv.eval.spline (., x, c (0, 1) )
	else
		.npuv.eval.direct (.npcdfuv.eval.ext, ., x)
}

npcdfuv.inverse = function (x, kernel.pdf=NA, kernel.cdf=sbc.cdf, nc=30, bw, smoothness=0.65, w=NA)
{  . = .npuv (TRUE, .npcdfuv.eval.ext, kernel.pdf, kernel.cdf, nc, smoothness, bw, w, x, TRUE)
	npcdfuv.f.inverse = function (y) {.npcdfuv.inverse.eval (y)}
  	attributes (npcdfuv.f.inverse) = list (
		class="npcdfuv.inverse",
		spline=TRUE, kernel.pdf=kernel.pdf, kernel.cdf=kernel.cdf, weighted=.$weighted,
		nc=.$nc, cy=.$cy, cx=.$cx, ct=.$ct, smoothness=.$smoothness, bw=.$bw, n=.$n, w=.$w, x=.$x)
    npcdfuv.f.inverse
}

.npcdfuv.inverse.eval = function (y)
{	. = attributes (sys.function (-1) )
	.npcdfuv.inverse.eval.spline (., y)
}

.npuv.eval.spline = function (., x, out)
{   n = length (x)
    y = numeric (n)
    for (i in 1:n)
    {   if (x [i] < .$cx [1])
			y [i] = out [1]
		else if (x [i] > .$cx [.$nc])
			y [i] = out [2]
        else
			y [i] = .spline.eval (.$nc, .$cx, .$cy, .$ct, x [i])
    }
    y
}

.npcdfuv.inverse.eval.spline = function (., y)
{  	n = length (y)
    x = numeric (n)
    for (i in 1:n)
    {   if (y [i] < 0 || y [i] > 1)
            stop ("y values must be between 0 and 1")
        else if (y [i] == 0)
            x [i] = .$cx [1]
        else if (y [i] == 1)
            x [i] = .$cx [.$nc]
		else
			x [i] = .spline.eval (.$nc, .$cy, .$cx, .$ct, y [i])
    }
    x
}

.npuv.eval.direct = function (ext, ., x)
{	n = length (x)
	y = numeric (n)
	for (i in 1:n)
		y [i] = ext (., x [i])
	y
}

.nppdfuv.eval.ext = function (., x)
{	y = 0
	for (i in 1:.$n)
	{	dist = x - .$x [i]
		y.local = 2 / .$bw * .$kernel.pdf (2 / .$bw * dist)
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

.npcdfuv.eval.ext = function (., x)
{	y = 0
	for (i in 1:.$n)
	{	dist = x - .$x [i]
		y.local = .$kernel.cdf (2 / .$bw * dist)
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
