npmean = function (nppdfuv.f)
{	. = attributes (nppdfuv.f)
	if (!inherits (nppdfuv.f, "nppdfuv") || !.$spline)
		stop ("needs nppdfuv object in spline form")
	y = 0
	for (i in 1:(.$nc - 1) )
	{	x = (.$cx [i] + .$cx [i + 1]) / 2
		y = y + .intt.integral.2 (.$cx [i], .$cx [i + 1], .$cy [i], .$cy [i + 1], .$ct [i], .$ct [i + 1], x)
	}
	y
}

npmode = function (nppdfuv.f, include.boundaries=TRUE, all=FALSE, warning=FALSE)
{	. = attributes (nppdfuv.f)
	if (!inherits (nppdfuv.f, "nppdfuv") || !.$spline)
		stop ("needs nppdfuv object in spline form")
	
	n = .$nc
	x = .$cx
	y = .$cy
	t = .$ct
	t.left = t [-n]
	t.right = t [-1]

	modes = values = numeric (0)
	if (include.boundaries)
	{	if (y [1] > y [2] && t.left [1] <=0)
		{	modes = x [1]
			values = y [1]
		}
		if (y [n - 1] < y [n] && t.right [n - 1] >=0)
		{	modes = c (modes, x [n])
			values = c (values, y [n])
		}
	}
	for (i in 1:(n - 1) )
	{	if ( (t.right [i] > 0 && t.left [i + 1] < 0) ||
			(t.right [i] >= 0 && t.left [i + 1] <= 0 && y [i] < y [i + 1] && y [i + 1] > y [i + 2]) ||
			(t.right [i] > 0 && t.left [i + 1] == 0 && y [i] < y [i + 1]) ||
			(t.right [i] == 0 && t.left [i + 1] < 0 && y [i] > y [i + 1]) )
		{	modes = c (modes, x [i + 1])
			values = c (values, t [i + 1])
		}
	}
	for (i in 1:(n - 1) )
	{	if ( (t.left [i] > 0 && t.right [i] < 0) ||
			(t.left [i] == 0 && t.right [i] < 0 && y [i] < y [i + 1]) ||
			(t.left [i] > 0 && t.right [i] == 0 && y [i] > y [i + 1]) )
		{	mode = .intt.argmax (x [i], x [i + 1], y [i], y [i + 1], t [i], t [i + 1])
			value = nppdfuv.f (mode)
			modes = c (modes, mode)
			values = c (values, value)
		}
	}
	mode.order = order (modes)
	modes = modes [mode.order]
	values = values [mode.order]

	if (warning)
	{	if (length (modes) == 0)
			warning ("no modal points")
		if (length (modes) > 1)
			warning ("multiple modal points")
	}
	if (all)
		modes
	else
		modes [which.max (values)]
}

nprng = function (npcdf.f.inverse, n)
{	. = attributes (npcdf.f.inverse)
	if (inherits (npcdf.f.inverse, "npcdfuv.inverse") )
		npcdf.f.inverse (runif (n) )
	else if (inherits (npcdf.f.inverse, "chained.npcdfmv.inverse") )
	{	x = matrix (runif (n * .$m), nrow=n)
		npcdf.f.inverse (x)
	}
	else
		stop ("\nerng needs\nnpcdfuv.inverse or chained.npcdfmv.inverse object")
}
