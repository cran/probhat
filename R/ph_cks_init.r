#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

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

.cksuv = function (f, classes, is.cdf, is.spline, nc, XLIM, bw.method, kernel, bw, smoothness, x, w, tail)
{   	xname = .varname (x)
	.test.nc (nc)
	x = .val.x.uv (x, TRUE)
	n = length (x)
	kernel = .val.k (kernel)
	is.weighted = (! (missing (w) || is.na (w [1]) ) )
	w = .val.w (is.weighted, n, w)
	xlim = range (x)
	if (missing (bw) )
		bw = auto.cbw (x, bw.method=bw.method, smoothness=smoothness)
	else
		smoothness = NA
	hbw = bw / 2
	xlim = xlim + c (-hbw, hbw)

	if (missing (tail) )
		.low = tail = NA
	else
	{	tail = .val.tail (tail)
		.low = (tail == "lower")
	}

	XLIM = .val.XLIM.uv (XLIM)
	.is.trunc = .is.trunc.uv (hbw, XLIM, xlim)
	.any.trunc = any (.is.trunc)

	.scalef = .xpnd = NA
	.any.trunc.lower = .is.trunc.lower = FALSE

	if (.any.trunc)
	{	.check.x.inside.uv (.is.trunc, XLIM, x)

		xlim = .update.xlim (.is.trunc, XLIM, xlim)
		.xpnd = .update.x.uv (.is.trunc, hbw, XLIM, n, x, w, is.weighted)
		if (is.cdf)
		{	.is.trunc.lower = .is.trunc.lower.side (rbind (.is.trunc), .low)
			.any.trunc.lower = any (.is.trunc.lower)
		}
		else
		{	subv = .pwith.eval.uv (bw, kernel@F, .xpnd$n, .xpnd$x, XLIM [1], XLIM [2],
				is.weighted, .xpnd$w)
			.scalef = 1 / subv
		}
	}
	data = list (n=n, xlim=xlim, x=x, w=w)

	sf = .EXTEND (f, classes,
		.const.cdf.lower=0,
		.any.trunc.lower, .is.trunc.lower,
		.any.trunc, .is.trunc, .xpnd, .scalef,
		.low,
		xname, XLIM, is.spline=FALSE, is.weighted, tail, spline.function=NA, kernel,
		bw, smoothness, n0=n, data)
	if (is.cdf)
	{	if (.any.trunc)
		{	if (.any.trunc.lower)
			{	lo = sf (.cdf.lower.side (.low, rbind (xlim) ), .ignore.trunc=TRUE)
				sf %$% ".const.cdf.lower" = lo
			}
			else
				lo = 0
			hi = sf (.cdf.upper.side (.low, rbind (xlim) ), .ignore.trunc=TRUE)
			sf %$% ".scalef" = 1 / (hi - lo)
		}
	}
	if (is.spline)
	{	sf %$% "spline.function" = .spline (sf, is.cdf, nc, .low)
		sf %$% "is.spline" = TRUE
	}
	sf
}

.cksmv = function (f, classes, XLIM, bw.method, kernel, bw, smoothness, x, w, init.class="mv", is.cdf=FALSE, tail)
{   	is.cond = (init.class == "c" || init.class == "mvc")
	xnames = .varnames (x, "x", is.cond)
	x = .val.x.mv (x)
	colnames (x) = xnames
	n = nrow (x)
	m = ncol (x)
	kernel = .val.k (kernel)
	is.weighted = (! (missing (w) || is.na (w [1]) ) )
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
	hbw = bw / 2
	xlim = xlim + outer (bw, c (-0.5, 0.5) )

	if (missing (tail) )
		.low = tail = NA
	else if (init.class == "mv")
	{	tail = .val.tail (tail, m)
		.low = (tail == "lower")
	}

	if (init.class != "ch")
	{	.any.trunc.lower = .is.trunc.lower = FALSE
		XLIM = .val.XLIM.mv (XLIM, m)
		.is.trunc = .is.trunc.mv (hbw, XLIM, xlim)
		.any.trunc = any (.is.trunc)
	
		.scalef = .xpnd = NA
		if (.any.trunc)
		{	.check.x.inside.mv (.is.trunc, XLIM, m, x)

			xlim = .update.xlim (.is.trunc, XLIM, xlim)
			.xpnd = .update.x.mv (.is.trunc, hbw, XLIM, n, m, x, w, is.weighted)

			if (init.class == "mv")
			{	if (is.cdf)
				{	.is.trunc.lower = .is.trunc.lower.side (.is.trunc, .low)
					.any.trunc.lower = any (.is.trunc.lower)
				}
				else
				{	subv = .pwith.eval.mv (bw, kernel@F, .xpnd$n, m, .xpnd$x, XLIM [,1], XLIM [,2],
						is.weighted, .xpnd$w)
					.scalef = 1 / subv
				}
			}
		}
	}
	data = list (n=n, xlim=xlim, x=x, w=w)

	sf = if (init.class == "mv")
		.EXTEND (f, classes,
			.any.trunc.lower, .is.trunc.lower,
			.any.trunc, .is.trunc, .scalef, .xpnd,
			.low,
			xnames, XLIM, is.weighted, tail, kernel, bw, smoothness, m, n0=n, data)
	else if (init.class == "c")
		.EXTEND (f, classes,
			.const.cdf.lower=0,
			.any.trunc.lower, .is.trunc.lower,
			.any.trunc, .is.trunc, .scalef, .xpnd,
			.constant=NA, .v=NA, .low=NA,
			M=NA, ncon=NA, xnames, XLIM,
			is.spline=FALSE, is.weighted, tail, spline.function=NA, conditions=NA,
			kernel, bw, smoothness, m, n0=0, data)
	else if (init.class == "mvc")
		.EXTEND (f, classes,
			.any.trunc.lower, .is.trunc.lower,
			.any.trunc, .is.trunc, .scalef, .xpnd,
			.constant=NA, .v=NA, .low=NA,
			M=NA, ncon=NA, xnames, XLIM,
			is.weighted, tail, conditions=NA,
			kernel, bw, smoothness, m, n0=0, data)
	else
		.EXTEND (f, classes,
			xnames, is.weighted, nc=NA, qf1=NA, kernel, bw, smoothness, m, n0=n, data)
	if (is.cdf && init.class == "mv")
	{	if (.any.trunc)
			sf %$% ".scalef" = .cdf.scalef (sf, m)
	}
	sf
}

.cksc = function (f, is.cdf, is.spline, nc, preserve.range, conditions, as.cset, warning, init.class="c")
{	m = f %$% "m"
	if (missing (conditions) )
		stop ("conditions required")
	conditions = rbind (conditions)
	ncon = ncol (conditions)
	if (ncon == 0)
		stop ("conditional models need at least one condition")
	M = m - ncon
	if (init.class == "c" && M != 1)
		stop ("\nuv-conditional models need one random variable\nor use mv-conditional models...")
	if (ncon < 0 || ncon > m - 1)
		stop ("length (conditions) not in interval [0, m - 1]")
	preserve.range = rep_len (preserve.range, M)
	if (is.cdf)
	{	tail = .val.tail (f %$% "tail", M)
		f %$% ".low" = (tail == "lower")
		f %$% "tail" = tail
	}
	if (ncon > 0)
	{	names = colnames (conditions)
		if (is.null (names) )
			names (conditions) = (f %$% "xnames")[1:ncon]
		else
		{	J = match (names, f %$% "xnames")
			if (any (is.na (J) ) )
				stop ("condition names not in variable names")
			J = c (J, (1:m)[-J])

			f %$% "xnames" = (f %$% "xnames")[J]
			f %$% "XLIM" = (f %$% "XLIM")[J,, drop=FALSE]
			f %$% "bw" = (f %$% "bw")[J]
			f %$% "smoothness" = (f %$% "smoothness")[J]

			if (f %$% ".any.trunc")
			{	f %$% ".is.trunc" = (f %$% ".is.trunc")[J,, drop=FALSE]
				f %$% ".xpnd" = .sub.data.col (f %$% ".xpnd", J)
			}
			f %$% "data" = .sub.data.col (f %$% "data", J)
		}
	}
	f %$% "M" = M
	f %$% "ncon" = ncon
	f %$% "conditions" = conditions
	f %$% "n0" = f %$% "n"

	CONS = 1:ncon
	RVS = (ncon + 1):m

	if (f %$% ".any.trunc")
	{	.check.x.inside.mv ( (f %$% ".is.trunc")[CONS,, drop=FALSE], (f %$% "XLIM")[CONS,, drop=FALSE], ncon, conditions, "condition")

		if (is.cdf)
		{	f %$% ".is.trunc.lower" = .is.trunc.lower.side ( (f %$% ".is.trunc")[RVS,, drop=FALSE], f %$% ".low")
			f %$% ".any.trunc.lower" = any (f %$% ".is.trunc.lower")
		}
	}

	nnull = 0
	nsets = nrow (conditions)

	VF = vector ("list", nsets)
	data0 = f %$% "data"
	for (ith.set in seq_len (nsets) )
	{	ith.conditions = conditions [ith.set,]
		data = f %$% "data"
		K = .con.sub (ncon, ith.conditions, f %$% "bw", data$n, data$x, data$w)
		if (K$n == 0)
			nnull = nnull + 1
		else
		{	new.f = f

			data$n = K$n
			data$x = K$x
			data$w = K$w
			xlim = data$xlim
			for (j in 1:m)
			{	if (j > ncon && preserve.range [j - ncon])
					0
				else
					xlim [j,] = range (data$x [,j]) + c (-0.5, 0.5) * (new.f %$% "bw")[j]
			}
			if (f %$% ".any.trunc")
				xlim = .update.xlim (f %$% ".is.trunc", f %$% "XLIM", xlim)
			data$xlim = xlim
			new.f %$% "data" = data

			if (f %$% ".any.trunc")
			{	xpnd = f %$% ".xpnd"
				K = .con.sub (ncon, ith.conditions, f %$% "bw", xpnd$n, xpnd$x, xpnd$w)
				xpnd$n = K$n
				xpnd$x = K$x
				xpnd$w = K$w
				new.f %$% ".xpnd" = xpnd
				fdata = xpnd

				. = attributes (new.f)
				if (is.cdf)
				{	if (init.class == "c")
					{	if (.$.low) u = .$XLIM [m, 1]
						else u = .$XLIM [m, 2]
						new.f %$% "const.cdf.lower" =
							.cdfuv.cks.eval.scalar (.$is.weighted, .$kernel@F, .$bw, xpnd$n, xpnd$x, xpnd$w, f %$% ".low", u)
					}
					else
						new.f %$% "scalef" = 0
				}
				else
				{	subv.top = .pwith.eval.mv (.$bw, .$kernel@F,
						xpnd$n, m, xpnd$x, .$XLIM [,1], .$XLIM [,2],
						.$is.weighted, xpnd$w)
					subv.btm = .pwith.eval.mv (.$bw [CONS], .$kernel@F,
						xpnd$n, ncon, xpnd$x [,CONS, drop=FALSE], .$XLIM [CONS, 1], .$XLIM [CONS, 2],
						.$is.weighted, xpnd$w)
					scalef = subv.btm / subv.top
					new.f %$% ".scalef" = scalef
				}
			}
			else
				fdata = data

			new.f %$% ".v" = .precompute.cksc.v (ncon, ith.conditions, (new.f %$% "kernel")@f, new.f %$% "bw", fdata$n, fdata$x)
			new.f %$% ".constant" = .sumk (new.f %$% "is.weighted", fdata$n, fdata$w, new.f %$% ".v")
			if (is.cdf)
			{	if (f %$% ".any.trunc")
				{	if (init.class == "c")
					{	low = f %$% ".low"
						if (f %$% ".any.trunc.lower")
						{	xlim = data$xlim
							lo = new.f (.cdf.lower.side (low, xlim [RVS,, drop=FALSE]), .ignore.trunc=TRUE)
						}
						else
							lo = 0
						hi = new.f (.cdf.upper.side (low, xlim [RVS,, drop=FALSE]), .ignore.trunc=TRUE)
						new.f %$% ".const.cdf.lower" = lo
						new.f %$% ".scalef" = 1 / (hi - lo)
					}
					else
						new.f %$% ".scalef" = .cdf.scalef (new.f, M)
				}
			}

			if (init.class == "c" && is.spline)
			{	new.f %$% "spline.function" = .spline (new.f, is.ccdf (new.f), nc, new.f %$% ".low")
				new.f %$% "is.spline" = TRUE
			}
			VF [[ith.set]] = new.f
		}
	}
	if (nsets == 1 && ! as.cset)
	{	if (nnull > 0 && warning)
			warning ("null model, no observations within conditional window")
		VF [[1]]
	}
	else
	{	if (nnull > 0 && warning)
			warning (sprintf ("%i null models, no observations within conditional windows", nnull) )
		VF
	}
}

.cksc.2 = function (f, classes, is.cdf,
	is.spline, nc, preserve.range, conditions, XLIM, bw.method, kernel, bw, smoothness, x, w, as.cset, warning, tail)
{	.test.nc (nc)
	f = .cksmv (f, classes, XLIM, bw.method, kernel, bw, smoothness, x, w, "c", FALSE, tail)
	.cksc (f, is.cdf, is.spline, nc, preserve.range, conditions, as.cset, warning)
}

.cksmvc.2 = function (f, classes, is.cdf,
	preserve.range, conditions, XLIM, bw.method, kernel, bw, smoothness, x, w, as.cset, warning, tail)
{	f = .cksmv (f, classes, XLIM, bw.method, kernel, bw, smoothness, x, w, "mvc", FALSE, tail)
	.cksc (f, is.cdf, NULL, NULL, preserve.range, conditions, as.cset, warning, "mvc")
}

.qfuv.cks = function (F, F.inv, nc, XLIM, bw.method, kernel, bw, smoothness, x, w, tail)
{	F = .cksuv (F, .CV.cdfuv.cks, TRUE, TRUE, nc, XLIM, bw.method, kernel, bw, smoothness, x, w, tail)
	attributes (F.inv) = attributes (F)
	F.inv %$% "class" = .CV.qfuv.cks
	F.inv %$% "spline.function" =
		.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy", F %$% ".low")
	F.inv
}

.qfc.cks = function (F, F.inv, nc, preserve.range, conditions, XLIM, bw.method, kernel, bw, smoothness, x, w, as.cset, warning, tail)
{	F = .cksc.2 (F, .CV.cdfc.cks, TRUE, TRUE, nc, preserve.range, conditions, XLIM, bw.method, kernel, bw, smoothness, x, w, as.cset, warning, tail)
	if (is.null (F) )
		NULL
	else if (is.ccdf (F) )
		.qfc.cks.ext (F, F.inv)
	else
	{	nsets = length (F)
		for (i in seq_len (nsets) )
		{	if (! is.null (F [[i]]) )
				F [[i]] = .qfc.cks.ext (F [[i]], F.inv)
		}
		F
	}
}

.qfc.cks.ext = function (F, F.inv)
{	attributes (F.inv) = attributes (F)
	F.inv %$% "class" = .CV.qfc.cks
	F.inv %$% "spline.function" =
		.modified.spline.transposed (F %$% "spline.function" %$% "cx", F %$% "spline.function" %$% "cy", F %$% ".low")
	F.inv
}

.chqf.cks = function (chqf.f, nc, bw.method, kernel, bw, smoothness, x, w)
{	.test.nc (nc)
	chqf.f = .cksmv (chqf.f, .CV.chqf.cks, NULL, bw.method, kernel, bw, smoothness, x, w, "ch")
	chqf.f %$% "nc" = nc
	chqf.f %$% "qf1" = .qfc.from.chqf.cks (chqf.f, 1)
	chqf.f
}

.qfc.from.chqf.cks = function (chqf.f, J, conditions=NA)
{	. = attributes (chqf.f)
	F = .cdfc4chqf.cks.eval
	F.inv = .qfc4chqf.cks.eval
	data = .$data
	xlim = data$xlim
	if (is.na (conditions [1]) )
		K = list (n = data$n, x = data$x, w = data$w)
	else
	{	K = .con.sub (J - 1, conditions, .$bw, data$n, data$x, data$w)
		xlim [J,] = range (K$x [,J]) + c (-0.5, 0.5) * .$bw [J]
	}
	data = list (
		xlim = xlim [1:J,, drop=FALSE],
		n = K$n,
		x = K$x [,1:J, drop=FALSE],
		w = K$w)
	attributes (F) = list (
		class = .CV.cdfc.cks,
		ncon = J - 1,
		is.weighted = .$is.weighted,
		nc = .$nc,
		kernel = .$kernel,
		bw = .$bw [1:J],
		m = J,
		data = data,
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
