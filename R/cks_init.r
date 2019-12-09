#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

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

.cksuv = function (f, superclass, subclass, is.spline, nc, kernel, bw, smoothness, x, w)
{   	variable.name = .varname (x)
	.test.nc (nc)
	x = .val.x.uv (x)
	n = length (x)
	kernel = .val.k (kernel)
	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w)
	xlim = range (x)
	if (missing (bw) )
		bw = smoothness * diff (xlim)
	else
		smoothness = bw / diff (xlim)
	xlim = xlim + c (-0.5, 0.5) * bw
	f = EXTEND (f, c (subclass, "cksuv", superclass, "cpd", "phmodel"),
		variable.name, is.spline=FALSE, is.weighted, spline.function=NA, kernel,
		bw, smoothness, xlim, n, x, w)
	if (is.spline)
	{	f %$% "spline.function" = .spline (f, superclass == "cdfuv", nc)
		f %$% "is.spline" = TRUE
	}
	f
}

.cksmv = function (f, classes, kernel, bw, smoothness, x, w, init.class="mv")
{   	variable.names = .varnames (x)
	x = .val.x.mv (x)
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
		for (j in 1:m)
			bw [j] = smoothness [j] * diff (xlim [j,])
	}
	else
	{	smoothness = numeric (m)
		bw = .val.params (m, bw)
		for (j in 1:m)
			smoothness [j] = bw [j] / diff (xlim [j,])
	}
	xlim = xlim + outer (bw, c (-0.5, 0.5) )
	if (init.class == "mv")
		EXTEND (f, classes,
			variable.names, is.weighted, kernel, bw, smoothness, xlim, n, m, x, w)
	else if (init.class == "c")
		EXTEND (f, classes,
			.constant=NA, .v=NA,
			M=NA, ncon=NA, variable.names,
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

.cksc = function (f, superclass.1, is.spline, nc, preserve.range, conditions, init.class="c")
{	m = f %$% "m"
	ncon = length (conditions)
	M = m - ncon
	if (init.class == "c" && M != 1)
		stop ("\nuv-conditional models need M = 1 and ncon = m - 1\nor use mv-conditional models...")
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
			f %$% "x" = (f %$% "x")[,J]
			f %$% "bw" = (f %$% "bw")[J]
			f %$% "smoothness" = (f %$% "smoothness")[J]
		}
	}
	f %$% "M" = M
	f %$% "ncon" = ncon
	f %$% "conditions" = conditions
	K = .con.sub (ncon, conditions, f %$% "bw", f %$% "n", f %$% "x", f %$% "w")
	if (K$n == 0)
		stop ("no observations within conditional region")
	f %$% "n" = K$n
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
	if (init.class == "c" && is.spline)
	{	f %$% "spline.function" = .spline (f, superclass.1 == "cdfuv", nc)
		f %$% "is.spline" = TRUE
	}
	f
}

.cksmv.2 = function (f, superclass, subclass, ...)
	.cksmv (f, c (subclass, "cksmv", superclass, "cpd", "phmodel"), ...)

.cksc.2 = function (f, superclass.1, superclass.2, subclass,
	is.spline, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
{	classes = c (subclass, "cksc", superclass.2, superclass.1, "cpd", "phmodel")
	.test.nc (nc)
	f = .cksmv (f, classes, kernel, bw, smoothness, x, w, "c")
	.cksc (f, superclass.1, is.spline, nc, preserve.range, conditions)
}

.cksmvc.2 = function (f, superclass.1, superclass.2, subclass,
	preserve.range, conditions, kernel, bw, smoothness, x, w)
{	classes = c (subclass, "cksmvc", superclass.2, superclass.1, "cpd", "phmodel")
	f = .cksmv (f, classes, kernel, bw, smoothness, x, w, "mvc")
	.cksc (f, NULL, NULL, NULL, preserve.range, conditions, "mvc")
}

.qfuv.cks = function (F, F.inv, nc, kernel, bw, smoothness, x, w)
{	F = .cksuv (F, "cdfuv", "cdfuv.cks", TRUE, nc, kernel, bw, smoothness, x, w)
	attributes (F.inv) = attributes (F)
	F.inv %$% "class" = c ("qfuv.cks", "cksuv", "qfuv", "cpd", "phmodel", "function")
	F.inv %$% "spline.function" =
		.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy")
	F.inv
}

.qfc.cks = function (F, F.inv, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
{	F = .cksc.2 (F, "cdfuv", "cdfc", "cdfc.cks", TRUE, nc, preserve.range, conditions, kernel, bw, smoothness, x, w)
	attributes (F.inv) = attributes (F)
	F.inv %$% "class" = c ("qfc.cks", "cksc", "qfc", "qfuv", "cpd", "phmodel", "function")
	F.inv %$% "spline.function" =
		.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy")
	F.inv 
}

.chqf.cks = function (chqf.f, nc, kernel, bw, smoothness, x, w)
{	.test.nc (nc)
	classes = c ("chqf.cks", "chcks", "chqf", "cpd", "phmodel")
	chqf.f = .cksmv (chqf.f, classes, kernel, bw, smoothness, x, w, "ch")
	chqf.f %$% "nc" = nc
	chqf.f %$% "qf1" = .qfc.from.chqf.cks (chqf.f, 1)
	chqf.f
}

.uv.from.mv.cks = function (f.mv, J)
{	. = attributes (f.mv)
	if (any (.$class == "pdfmv") )
	{	f = .pdfuv.cks.eval
		class = c ("pdfuv.cks", "pdfuv")
	}
	else
	{	f = .cdfuv.cks.eval
		class = c ("cdfuv.cks", "cdfuv")
	}
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
	f %$% "spline.function" = .spline (f, class [2] == "cdfuv", 30)
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
		class = c ("cdfc.cks", "cdfuv", "phmodel"),
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
	class (F.inv) = c ("qfc.cks", "qfuv", "phmodel")
	cx = seq (F, .$nc)
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
