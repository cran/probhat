#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.test.nc = function (nc)
{	if (nc < 10)
		stop ("spline models need >= 10 control points")
}

.cksuv = function (f, classes, is.cdf, is.spline, nc, bw.method, kernel, bw, smoothness, x, w)
{   	variable.name = .varname (x)
	.test.nc (nc)
	x = .val.x.uv (x, TRUE)
	n = length (x)
	kernel = .val.k (kernel)
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w)
	xlim = range (x)
	if (missing (bw) )
		bw = auto.cbw (x, bw.method=bw.method, smoothness=smoothness)
	else
		smoothness = NA
	xlim = xlim + c (-0.5, 0.5) * bw
	f = EXTEND (f, classes,
		variable.name, is.spline=FALSE, is.weighted, spline.function=NA, kernel,
		bw, smoothness, xlim, n, x, w)
	if (is.spline)
	{	f %$% "spline.function" = .spline (f, is.cdf, nc)
		f %$% "is.spline" = TRUE
	}
	f
}

.cksmv = function (f, classes, bw.method, kernel, bw, smoothness, x, w, init.class="mv")
{   	variable.names = .varnames (x)
	x = .val.x.mv (x)
	colnames (x) = variable.names
	n = nrow (x)
	m = ncol (x)
	kernel = .val.k (kernel)
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w)
	xlim = matrix (0, m, 2)
	for (j in 1:m)
		xlim [j,] = range (x [,j])
	if (missing (bw) )
	{	bw = numeric (m)
		smoothness = .val.params (m, smoothness)
		bw = auto.cbw (x, bw.method=bw.method, smoothness=smoothness)
	}
	else
	{	bw = .val.params (m, bw)
		smoothness = rep (NA, m)
	}
	xlim = xlim + outer (bw, c (-0.5, 0.5) )
	if (init.class == "mv")
		EXTEND (f, classes,
			variable.names, is.weighted, kernel, bw, smoothness, xlim, n, m, x, w)
	else if (init.class == "c")
		EXTEND (f, classes,
			.constant=NA, .v=NA,
			M=NA, ncon=NA, variable.name=NA, variable.names,
			is.spline=FALSE, is.weighted, spline.function=NA, conditions=NA,
			kernel, bw, smoothness, xlim, n, m, x, w)
	else if (init.class == "mvc")
		EXTEND (f, classes,
			.constant=NA, .v=NA,
			M=NA, ncon=NA, variable.names,
			is.weighted, conditions=NA,
			kernel, bw, smoothness, xlim, n, m, x, w)
	else
		EXTEND (f, classes,
			variable.names, is.weighted, nc=NA, qf1=NA, kernel, bw, smoothness, xlim, n, m, x, w)
}

.cksc = function (f, is.cdf, is.spline, nc, preserve.range, conditions, warning, init.class="c")
{	m = f %$% "m"
	if (missing (conditions) )
		stop ("conditions required")
	ncon = length (conditions)
	if (ncon == 0)
		stop ("conditional models need at least one condition")
	M = m - ncon
	if (init.class == "c" && M != 1)
		stop ("\nuv-conditional models need one random variable\nor use mv-conditional models...")
	if (ncon < 0 || ncon > m - 1)
		stop ("length (conditions) not in interval [0, m - 1]")
	else if (ncon > 0)
	{	names = names (conditions)
		if (is.null (names) )
			names (conditions) = (f %$% "variable.names")[1:ncon]
		else
		{	J = match (names, f %$% "variable.names")
			if (any (is.na (J) ) )
				stop ("condition names not in variable names")
			J = c (J, (1:m)[-J])
			f %$% "variable.names" = (f %$% "variable.names")[J]
			f %$% "xlim" = (f %$% "xlim")[J,]
			f %$% "bw" = (f %$% "bw")[J]
			f %$% "smoothness" = (f %$% "smoothness")[J]
			f %$% "x" = (f %$% "x")[,J]
		}
	}
	f %$% "M" = M
	f %$% "ncon" = ncon
	f %$% "conditions" = conditions
	K = .con.sub (ncon, conditions, f %$% "bw", f %$% "n", f %$% "x", f %$% "w")
	if (K$n == 0)
	{	warning ("no observations within conditional window")
		NULL
	}
	else
	{	f %$% "n" = K$n
		f %$% "x" = K$x
		f %$% "w" = K$w
		if (! preserve.range)
		{	if (init.class == "c")
				xlim = range (K$x [,m]) + c (-0.5, 0.5) * (f %$% "bw")[m]
			else
			{	xlim = matrix (0, M, 2)
				for (j in 1:M)
					xlim [j,] = range ( (f %$% "x")[,ncon + j]) + c (-0.5, 0.5) * (f %$% "bw")[ncon + j]
			}
			f %$% "xlim" = xlim
		}
		else
			f %$% "xlim" = (f %$% "xlim")[(ncon + 1):m,]
		f %$% ".v" = .precompute.cksc.v (ncon, conditions, (f %$% "kernel")$pdf, f %$% "bw", f %$% "n", f %$% "x")
		f %$% ".constant" = .sumk (f %$% "is.weighted", f %$% "n", f %$% "w", f %$% ".v")
		if (init.class == "c")
		{	f %$% "variable.name" = (f %$% "variable.names")[m]
			f %$% "variable.names" = NULL
			if (is.spline)
			{	f %$% "spline.function" = .spline (f, is.ccdf (f), nc)
				f %$% "is.spline" = TRUE
			}
		}
		f
	}
}

.cksc.2 = function (f, classes, is.cdf,
	is.spline, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
{	.test.nc (nc)
	f = .cksmv (f, classes, bw.method, kernel, bw, smoothness, x, w, "c")
	.cksc (f, is.cdf, is.spline, nc, preserve.range, conditions, warning)
}

.cksmvc.2 = function (f, classes,
	preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
{	f = .cksmv (f, classes, bw.method, kernel, bw, smoothness, x, w, "mvc")
	.cksc (f, NULL, NULL, NULL, preserve.range, conditions, warning, "mvc")
}

.qfuv.cks = function (F, F.inv, nc, bw.method, kernel, bw, smoothness, x, w)
{	F = .cksuv (F, .CV.cdfuv.cks, TRUE, TRUE, nc, bw.method, kernel, bw, smoothness, x, w)
	attributes (F.inv) = attributes (F)
	F.inv %$% "class" = .CV.qfuv.cks
	F.inv %$% "spline.function" =
		.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy")
	F.inv
}

.qfc.cks = function (F, F.inv, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
{	F = .cksc.2 (F, .CV.cdfuv.cks, TRUE, TRUE, nc, preserve.range, conditions, bw.method, kernel, bw, smoothness, x, w, warning)
	if (is.null (F) )
		F
	else
	{	attributes (F.inv) = attributes (F)
		F.inv %$% "class" = .CV.qfc.cks
		F.inv %$% "spline.function" =
			.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy")
		F.inv
	}
}

.chqf.cks = function (chqf.f, nc, bw.method, kernel, bw, smoothness, x, w)
{	.test.nc (nc)
	chqf.f = .cksmv (chqf.f, .CV.chqf.cks, bw.method, kernel, bw, smoothness, x, w, "ch")
	chqf.f %$% "nc" = nc
	chqf.f %$% "qf1" = .qfc.from.chqf.cks (chqf.f, 1)
	chqf.f
}

.uv.from.mv.cks = function (f.mv, J)
{	if (is.pdfmv (f.mv) && is.cks (f.mv) )
	{	f = .pdfuv.cks.eval
		class = c ("PRIVATE_OBJECT", .CV.pdfuv.cks)
	}
	else if (is.ccdfmv (f.mv) && is.cks (f.mv) )
	{	f = .cdfuv.cks.eval
		class = c ("PRIVATE_OBJECT", .CV.cdfuv.cks)
	}
	else
		stop ("error in .uv.from.mv.cks")

	. = attributes (f.mv)
	attributes (f) = list (
		class = class,
		variable.name = .$variable.names [J],
		is.spline = FALSE,
		is.weighted = .$is.weighted,
		spline.function = NA,
		kernel = .$kernel,
		bw = .$bw [J],
		xlim = .$xlim [J,],
		n = .$n,
		x = .$x [,J],
		w = .$w)
	f %$% "spline.function" = .spline (f, is.ccdf (f), 30)
	f %$% "is.spline" = TRUE
	f
}

.qfc.from.chqf.cks = function (chqf.f, J, conditions=NA)
{	. = attributes (chqf.f)
	F = .cdfc.cks.eval.2
	F.inv = .qfc.cks.eval.2
	if (is.na (conditions [1]) )
	{	K = list (n = .$n, x = .$x, w = .$w)
		xlim = .$xlim [J,]
	}
	else
	{	K = .con.sub (J - 1, conditions, .$bw, .$n, .$x, .$w)
		xlim = range (K$x [,J]) + c (-0.5, 0.5) * .$bw [J]
	}
	attributes (F) = list (
		class = .CV.cdfc.cks,
		ncon = J - 1,
		is.weighted = .$is.weighted,
		nc = .$nc,
		kernel = .$kernel,
		bw = .$bw [1:J],
		xlim = xlim,
		n = K$n,
		x = K$x [,1:J, drop=FALSE],
		w = K$w,
		conditions = conditions)
	attributes (F.inv) = attributes (F)
	class (F.inv) = .CV.qfc.cks
	cx = seq (F, n = .$nc)
	F.inv %$% "spline.function" = .modified.spline.transposed (cx, F (cx) )
	F.inv
}

.con.sub = function (ncon, conditions, bw, n, x, w)
{	I = rep (TRUE, n)
	for (j in seq_len (ncon) )
	{	xlim = conditions [j] + c (-0.5, 0.5) * bw [j]
		I = (I & x [,j] >= xlim [1] & x [,j] <= xlim [2])
	}
	if (! is.na (w [1]) )
		w = w [I]
	list (n = sum (I), x = x [I,, drop=FALSE], w=w)
}
