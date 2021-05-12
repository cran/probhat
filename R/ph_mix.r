#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

#unpacked by .gmix, .xmix
.mix = function (continuous, g, x, w, conditions, throw.warning, ...)
{	objs = .cat.data (g, TRUE)
	gnames = nlevels = levels = n = m = g = NULL
	.UNPACK (objs)

	mg = m
	xnames = .varnames (x, "x", TRUE)
	x = .val.x.uv.or.mv (x)
	colnames (x) = xnames
	mx = ncol (x)

	is.weighted = (! (missing (w) || is.na (w [1]) ) )
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
	{	cstr = paste (connames, collapse=" ")
		bstr = paste (bnames, collapse=" ")
		err = paste0 (
			"condition names not in variable names:\n    ",
				"cons: ", cstr, "\n    ",
				"vars: ", bstr)
		stop (err)
	}

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

	n0 = n
	if (ngcon > 0)
	{	gJ = c (gJ, (1:mg)[-gJ])
		gnames = gnames [gJ]
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
		xnames = xnames [xJ]
		x = x [,xJ, drop=FALSE]
	}

	if (n == 0)
	{	if (throw.warning)
			warning ("no observations within conditional window")
		NULL
	}
	else
	{	if (continuous)
			.LIST (ngcon, nxcon, gnames, xnames, gconditions, xconditions, n0, mg, mx, x, w)
		else
		{	variable.name = gnames [mg]
			nlevels = nlevels [mg]
			levels = levels [[mg]]
			g = g [,mg]
			.LIST (is.weighted, gnames, xnames, nlevels, levels, n0, n, mg, mx, g, x, w, gconditions, xconditions)
		}
	}
}

.gmix = function (f, cv, g, x, w, conditions, throw.warning, ..., bw, smoothness=1)
{	mixobj = .mix (FALSE, g, x, w, conditions, throw.warning, ...)
	if (is.null (mixobj) )
		NULL
	else
	{	is.weighted = gnames = xnames = nlevels = levels = n0 = n = mg = mx = g = x = w = gconditions = xconditions = NULL
		.UNPACK (mixobj)

		if (is.weighted)
			h = w
		else
			h = rep (1, n)
		.gsum = sum (h)
	
		.probs = .iterate.uv (.pmfuv.cat.eval.scalar, n, g, h / sum (h), u=1:nlevels)
		.PROBS = cumsum (.probs)
		.PROBS [nlevels] = 1
		gfh = .EXTEND (.pmfuv.cat.eval, .CV.pmfuv.cat,
			.probs, .PROBS,
			freq=FALSE, nlevels, glim = c (1, nlevels), levels, n, g, h)

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
			.EXTEND (f, cv,
				.probs, .PROBS, .gsum,
				gnames, xnames,
				gconditions, xconditions, nlevels, glim = c (1, nlevels), levels, n, n0, mg, mx, g, h)
		}
	}
}

.xmix = function (constructor, constructor.c, cv, g, x, w, conditions, throw.warning, ...)
{	mixobj = .mix (TRUE, g, x, w, conditions, throw.warning)
	if (is.null (mixobj) )
		NULL
	else
	{	ngcon = nxcon = NULL
		gnames = xnames = gconditions = xconditions = n0 = mg = mx = x = w = NULL
		.UNPACK (mixobj)

		if (length (xconditions) == 0)
		{	sf = constructor (x, ..., w=w)
		}
		else
			sf = constructor.c (x, ..., conditions=xconditions, w=w)
		class (sf) = cv
		sf %$% ".ngcon" = ngcon
		sf %$% ".nxcon" = nxcon
		sf %$% "gnames" = gnames
		sf %$% "xnames" = xnames
		sf %$% "gconditions" = gconditions
		sf %$% "xconditions" = xconditions
		sf %$% "n0" = n0
		sf %$% "mg" = mg
		sf %$% "mx" = mx
		sf
	}
}

ph4.pmfc.gmix = function (g, x, ..., conditions, warning=TRUE, w)
{	.arg.error (...)
	.gmix (.pmfuv.cat.eval, .CV.pmfc.gmix, g, x, w, conditions, warning, ...)
}

ph4.cdfc.gmix = function (g, x, ..., conditions, warning=TRUE, w)
{	.arg.error (...)
	.gmix (.cdfuv.cat.eval, .CV.cdfc.gmix, g, x, w, conditions, warning, ...)
}

ph4.qfc.gmix = function (g, x, ..., conditions, warning=TRUE, w)
{	.arg.error (...)
	.gmix (.qfuv.cat.eval, .CV.qfc.gmix, g, x, w, conditions, warning, ...)
}

ph4.pdfc.xmix = function (g, x, ..., conditions, warning=TRUE, w)
	.xmix (pdfuv.cks, pdfc.cks, .CV.pdfc.xmix, g, x, w, conditions, warning, ...)

ph4.cdfc.xmix = function (g, x, ..., conditions, warning=TRUE, w)
	.xmix (cdfuv.cks, cdfc.cks, .CV.cdfc.xmix, g, x, w, conditions, warning, ...)

ph4.qfc.xmix = function (g, x, ..., conditions, warning=TRUE, w)
	.xmix (qfuv.cks, qfc.cks, .CV.qfc.xmix, g, x, w, conditions, warning, ...)
