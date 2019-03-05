.eccdf = function (w, x, inverse=FALSE)
{   x = as.numeric (x)
	n = length (x)
	if (any (is.na (x) ) )
        stop ("no missing values allowed")
	if (length (unique (x) ) < 3)
        stop ("x needs 3 or more unique values")

	weighted = (!is.na (w [1]) )
	if (weighted)
	{	w = as.numeric (w)
		if (n != length (w) )
			stop ("length(w) not equal length(x)")
		if (any (is.na (w) ) )
			stop ("no missing values allowed")
		if (round (sum (w), 1) != 1)
			stop ("sum(w) must be approx 1")
		if (any (w <= 0) )
			stop ("all w need to be > 0")
	}

    all.unique = (n == length (unique (x) ) )
    if (!all.unique)
    {   sd = sd (x) / 1000
		while (!all.unique)
		{   x = x + runif (n, -sd, sd)
			all.unique = (n == length (unique (x) ) )
		}
    }

	x.order = order (x)
	x = x [x.order]

	if (weighted)
	{	w = w [x.order]
		intw = (w [-n] + w [-1]) / 2
		outw = (w [1] + w [n]) / 2
		intw =  intw / (1 - outw)
		y = cumsum (c (0, intw) )
	}
	else
		y = ( (1:n) - 1) / (n - 1)
	y [1] = 0
	y [n] = 1
	if (inverse)
		t = .slopes (n, y, x)
	else
		t = .slopes (n, x, y)

	list (n=n, w=w, x=x, y=y, t=t)
 }
 
eccdf = function (x, w=NA)
{   . = .eccdf (w, x)
	eccdf.f = function (x) {.eccdf.eval (x)}
	attributes (eccdf.f) = list (class="eccdf", n=.$n, w=.$w, x=.$x, y=.$y, t=.$t)
	eccdf.f
}

.eccdf.eval = function (x)
{   . = attributes (sys.function (-1) )
	n = length (x)
    y = numeric (n)
    for (i in 1:n)
    {   if (x [i] < .$x [1])
			y [i] = 0
		else if (x [i] > .$x [.$n])
			y [i] = 1
        else
			y [i] = .spline.eval (.$n, .$x, .$y, .$t, x [i])
    }
    y
}

eccdf.inverse = function (x, w=NA)
{	. = .eccdf (w, x, TRUE)
	eccdf.f.inverse = function (y) {.eccdf.inverse.eval (y)}
  	attributes (eccdf.f.inverse) = list (class="eccdf.inverse", n=.$n, w=.$w, y=.$y, x=.$x, t=.$t)
    eccdf.f.inverse
}

.eccdf.inverse.eval = function (y)
{   . = attributes (sys.function (-1) )
    n = length (y)
    x = numeric (n)
    for (i in 1:n)
    {   if (y [i] < 0 || y [i] > 1)
            stop ("y values must be between 0 and 1")
		x [i] = .spline.eval (.$n, .$y, .$x, .$t, y [i])
    }
    x
}
