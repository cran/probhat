#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

#unpacked by .gmix, .xmix
.mix = function (continuous, g, x, w, conditions, throw.warning, ...)
{	objs = .cat.data (g)
	variable.names = nlevels = levels = n = m = g = NULL
	UNPACK (objs)

	gnames = variable.names
	mg = m

	xnames = .varnames (x)
	x = .val.x.uv.or.mv (x)
	colnames (x) = xnames
	mx = ncol (x)

	is.weighted = (! is.na (w [1]) )
	w = .val.w (is.weighted, n, w, FALSE)

	if (mg == 0 || mx == 0)
		stop ("need at least one categorical and continuous variable")
	if (n != nrow (x) )
		stop ("number of categorical obs != number continuous obs")
	if (any (gnames %in% xnames) )
		stop ("categorical and continuous variables need unique names")

	ncon = length (conditions)
	bnames = c (gnames, xnames)
	connames = names (conditions)
	if (is.null (connames) )
		stop ("unnamed conditions in mixed conditional model")
	if (! all (connames %in% bnames) )
		stop ("condition names not in variable names")

	gconditions = conditions [connames %in% gnames]
	xconditions = unlist (conditions [connames %in% xnames])
	ngcon = length (gconditions)
	nxcon = length (xconditions)

	if (continuous)
	{	if (ngcon != mg)
			stop ("number of g conditions != number of g variables")
		if (nxcon != mx - 1)
			stop ("number of x conditions != number of x variables minus one")
	}
	else
	{	if (ngcon != mg - 1)
			stop ("number of g conditions != number of g variables minus one")
		if (nxcon != mx)
			stop ("number of x conditions != number of x variables")
	}

	gJ = match (names (gconditions), gnames)
	xJ = match (names (xconditions), xnames)
	if (ngcon > 0)
	{	gJ = c (gJ, (1:mg)[-gJ])
		variable.names = variable.names [gJ]
		nlevels = nlevels [gJ]
		levels = levels [gJ]
		g = g [,gJ, drop=FALSE]

		I = rep (TRUE, n)
		for (k in seq_len (ngcon) )
			I = I & .in.col (nlevels [k], levels [[k]], gconditions [[k]], g [,k])
		n = sum (I)

		g = g [I,, drop=FALSE]
		x = x [I,, drop=FALSE]
		w = w [I]
	}
	if (nxcon > 0)
	{	xJ = c (xJ, (1:mx)[-xJ])
		x = x [,xJ, drop=FALSE]
	}

	if (n == 0)
	{	if (throw.warning)
			warning ("no observations within conditional window")
		NULL
	}
	else
	{	if (continuous)
			LIST (x, w, gconditions, xconditions)
		else
		{	variable.name = gnames [mg]
			nlevels = nlevels [mg]
			levels = levels [[mg]]
			g = g [,mg]
			LIST (is.weighted, variable.name, nlevels, levels, n, g, x, w, gconditions, xconditions)
		}
	}
}

.gmix = function (f, cv, g, x, w, conditions, throw.warning, ..., bw, smoothness=1, freq=FALSE)
{	mixobj = .mix (FALSE, g, x, w, conditions, throw.warning, ...)
	if (is.null (mixobj) )
		NULL
	else
	{	is.weighted = variable.name = nlevels = levels = n = g = x = w = gconditions = xconditions = NULL
		UNPACK (mixobj)
	
		.probs = .iterate.uv (.pmfuv.cat.eval.scalar, is.weighted, n, g, w / sum (w), u=1:nlevels)
		.PROBS = cumsum (.probs)
		.PROBS [nlevels] = 1
		gfh = EXTEND (.pmfuv.cat.eval, .CV.pmfuv.cat,
			.probs, .PROBS,
			variable.name, is.weighted, freq=FALSE, nlevels, xlim = c (1, nlevels), levels, n, g, w)

		bw = auto.cbw (x, bw=bw, smoothness=smoothness)
		cfh = pdfmv.cks (x, ..., bw=bw, w=w)
		const = cfh (xconditions)

		if (const <= 0)
		{	if (throw.warning)
				warning ("no observations within conditional window")
			NULL
		}
		else
		{	cfh.set = vector ("list", nlevels)
			for (k in 1:nlevels)
			{	I = (g == k)
				if (sum (I) > 0)
				{	xsub = x [I,, drop=FALSE]
					wsub = NA
					if (is.weighted)
						wsub = w [I]
					cfh.set [[k]] = pdfmv.cks (xsub, ..., bw=bw, w=wsub)
				}
			}
	
			.probs = rep (nlevels, 0)
			for (k in 1:nlevels)
			{	if (! is.null (cfh.set [[k]]) )
					.probs [k] = gfh (k) * cfh.set [[k]](xconditions) / const
			}
	
			.PROBS = cumsum (.probs)
			.PROBS [nlevels] = 1
			EXTEND (f, cv,
				.probs, .PROBS,
				is.weighted, variable.name, gconditions, xconditions, freq, nlevels, xlim = c (1, nlevels), levels, n, g, w)
		}
	}
}

.xmix = function (constructor, constructor.c, g, x, w, conditions, throw.warning, ...)
{	mixobj = .mix (TRUE, g, x, w, conditions, throw.warning)
	if (is.null (mixobj) )
		NULL
	else
	{	x = w = gconditions = xconditions = NULL
		UNPACK (mixobj)

		if (length (xconditions) == 0)
		{	sf = constructor (x, ..., conditions=xconditions, w=w)
			sf = EXTEND (sf, "cpdc")
		}
		else
			sf = constructor (x, ..., conditions=xconditions, w=w)
		sf %$% "gconditions" = gconditions
		sf %$% "xconditions" = xconditions
		sf
	}
}

pmfc.gmixp = function (g, x, ..., conditions, warning=TRUE, w=NA, freq=FALSE)
	.gmix (.pmfuv.cat.eval, .CV.pmfc.cat, g, x, w, conditions, warning, ..., freq=freq)

cdfc.gmixp = function (g, x, ..., conditions, warning=TRUE, w=NA)
	.gmix (.cdfuv.cat.eval, .CV.cdfc.cat, g, x, w, conditions, warning, ...)

qfc.gmixp = function (g, x, ..., conditions, warning=TRUE, w=NA)
	.gmix (.qfuv.cat.eval, .CV.qfc.cat, g, x, w, conditions, warning, ...)

pdfc.xmixp = function (g, x, ..., conditions, warning=TRUE, w=NA)
	.xmix (pdfuv.cks, pdfc.cks, g, x, w, conditions, warning, ...)

cdfc.xmixp = function (g, x, ..., conditions, warning=TRUE, w=NA)
	.xmix (cdfuv.cks, cdfc.cks, g, x, w, conditions, warning, ...)

qfc.xmixp = function (g, x, ..., conditions, warning=TRUE, w=NA)
	.xmix (qfuv.cks, qfc.cks, g, x, w, conditions, warning, ...)
