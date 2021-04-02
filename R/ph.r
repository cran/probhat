#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.CV = c ("phmodel", "phpd", "phob")

.CV.pmfuv.dks = c ("pmfuv.dks", "dksuv", "dks", "pmfuv", "dpduv", "pmf", "dpd", .CV)
.CV.cdfuv.dks = c ("cdfuv.dks", "dksuv", "dks", "dcdfuv", "dpduv", "dcdf", "dpd", .CV)
.CV.qfuv.dks = c ("qfuv.dks", "dksuv", "dks", "dqfuv", "dpduv", "dqf", "dpd", .CV)

.CV.pdfuv.cks = c ("pdfuv.cks", "cksuv", "cks", "pdfuv", "cpduv", "pdf", "cpd", .CV)
.CV.cdfuv.cks = c ("cdfuv.cks", "cksuv", "cks", "ccdfuv", "cpduv", "ccdf", "cpd", .CV)
.CV.qfuv.cks = c ("qfuv.cks", "cksuv", "cks", "cqfuv", "cpduv", "cqf", "cpd", .CV)
.CV.pdfmv.cks = c ("pdfmv.cks", "cksmv", "cks", "pdfmv", "cpdmv", "pdf", "cpd", .CV)
.CV.cdfmv.cks = c ("cdfmv.cks", "cksmv", "cks", "ccdfmv", "cpdmv", "ccdf", "cpd", .CV)
.CV.pdfc.cks = c ("pdfc.cks", "cksc", "cks", "pdfc", "cpdc", "pdfuv", "cpduv", "pdf", "cpd", .CV)
.CV.cdfc.cks = c ("cdfc.cks", "cksc", "cks", "ccdfc", "cpdc", "ccdfuv", "cpduv", "ccdf", "cpd", .CV)
.CV.qfc.cks = c ("qfc.cks", "cksc", "cks", "cqfc", "cpdc", "cqfuv", "cpduv", "cqf", "cpd", .CV)
.CV.pdfmvc.cks = c ("pdfmvc.cks", "cksmvc", "cks", "pdfmvc", "cpdmvc", "pdfmv", "cpdmv", "pdf", "cpd", .CV)
.CV.cdfmvc.cks = c ("cdfmvc.cks", "cksmvc", "cks", "ccdfmvc", "cpdmvc", "ccdfmv", "cpdmv", "ccdf", "cpd", .CV)
.CV.chqf.cks = c ("cchqf.cks", "cks", "cchqf", "cpd", .CV)

.CV.pmfuv.cat = c ("pmfuv.cat", "catuv", "cat", "pmfuv", "dpduv", "pmf", "dpd", .CV)
.CV.cdfuv.cat = c ("cdfuv.cat", "catuv", "cat", "dcdfuv", "dpduv", "dcdf", "dpd", .CV)
.CV.qfuv.cat = c ("qfuv.cat", "catuv", "cat", "dqfuv", "dpduv", "dqf", "dpd", .CV)
.CV.pmfc.cat = c ("pmfc.cat", "catc", "cat", "pmfc", "pmfuv", "dpdc", "dpduv", "pmf", "dpd", .CV)
.CV.cdfc.cat = c ("cdfc.cat", "catc", "cat", "dcdfc", "dcdfuv", "dpdc", "dpduv", "dcdf", "dpd", .CV)
.CV.qfc.cat = c ("qfc.cat", "catc", "cat", "dqfc", "dqfuv", "dpdc", "dpduv", "dqf", "dpd", .CV)

.CV.cdfuv.el = c ("cdfuv.el", "eluv", "el", "ccdfuv", "cpduv", "ccdf", "cpd", .CV)
.CV.qfuv.el = c ("qfuv.el", "eluv", "el", "cqfuv", "cpduv", "cqf", "cpd", .CV)

.CV.pmfc.gmix = c ("pmfc.gmix", "gmixc", "gmix", "pmfc", "pmfuv", "dpdc", "dpduv", "pmf", "dpd", .CV)
.CV.cdfc.gmix = c ("cdfc.gmix", "gmixc", "gmix", "dcdfc", "dcdfuv", "dpdc", "dpduv", "dcdf", "dpd", .CV)
.CV.qfc.gmix = c ("qfc.gmix", "gmixc", "gmix", "dqfc", "dqfuv", "dpdc", "dpduv", "dqf", "dpd", .CV)

.CV.pdfc.xmix = c ("pdfc.xmix", "xmixc", "xmix", "pdfc", "cpdc", "pdfuv", "cpduv", "pdf", "cpd", .CV)
.CV.cdfc.xmix = c ("cdfc.xmix", "xmixc", "xmix", "ccdfc", "cpdc", "ccdfuv", "cpduv", "ccdf", "cpd", .CV)
.CV.qfc.xmix = c ("qfc.xmix", "xmixc", "xmix", "cqfc", "cpdc", "cqfuv", "cpduv", "cqf", "cpd", .CV)

.inc = function (sf, include, class0, class1)
{	if (include) (inherits (sf, class0) || inherits (sf, class1) )
	else (inherits (sf, class0) )
}

pmf.dks = function (...) pmfuv.dks (...)
cdf.dks = function (...) cdfuv.dks (...)
qf.dks = function (...) qfuv.dks (...)
pdf.cks = function (...) pdfuv.cks (...)
cdf.cks = function (...) cdfuv.cks (...)
qf.cks = function (...) qfuv.cks (...)
pmf.cat = function (...) pmfuv.cat (...)
cdf.cat = function (...) cdfuv.cat (...)
qf.cat = function (...) qfuv.cat (...)
cdf.el  = function (...) cdfuv.el (...)
qf.el = function (...) qfuv.el (...)

is.phob = function (object) inherits (object, "phob")
is.phpd = function (object) inherits (object, "phpd")
is.phmodel = function (object) inherits (object, "phmodel")

is.dpd = function (object) inherits (object, "dpd")
is.cpd = function (object) inherits (object, "cpd")

is.pmf = function (object) inherits (object, "pmf")
is.dcdf = function (object) inherits (object, "dcdf")
is.dqf = function (object) inherits (object, "dqf")
is.pdf = function (object) inherits (object, "pdf")
is.ccdf = function (object) inherits (object, "ccdf")
is.cqf = function (object) inherits (object, "cqf")
is.cchqf = function (object) inherits (object, "cchqf")

is.dpduv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "dpduv", "dpdc")
is.dpdc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "dpdc", "dpdmvc")
is.cpduv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "cpduv", "cpdc")
is.cpdmv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "cpdmv", "cpdmvc")
is.cpdc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "cpdc", "cpdmvc")
is.cpdmvc = function (object) inherits (object, "cpdmvc")

is.pduv = function (object, include.conditional=TRUE)
	(is.dpduv (object, include.conditional) || is.cpduv (object, include.conditional) )

is.pmfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "pmfuv", "pmfc")
is.pmfc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "pmfc", "pmfmvc")
is.dcdfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "dcdfuv", "dcdfc")
is.dcdfc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "dcdfc", "dcdfmvc")
is.dqfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "dqfuv", "dqfc")
is.dqfc = function (object) inherits (object, "dqfc")
is.pdfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "pdfuv", "pdfc")
is.pdfmv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "pdfmv", "pdfmvc")
is.pdfc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "pdfc", "pdfmvc")
is.pdfmvc = function (object) inherits (object, "pdfmvc")
is.ccdfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "ccdfuv", "ccdfc")
is.ccdfmv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "ccdfmv", "ccdfmvc")
is.ccdfc = function (object, include.multivariate=TRUE) .inc (object, include.multivariate, "ccdfc", "ccdfmvc")
is.ccdfmvc = function (object) inherits (object, "ccdfmvc")
is.cqfuv = function (object, include.conditional=TRUE) .inc (object, include.conditional, "cqfuv", "cqfc")
is.cqfc = function (object) inherits (object, "cqfc")

is.phspline = function (object)
{	if (is.pdfuv (object, TRUE) || is.ccdfuv (object, TRUE) )
		object %$% "is.spline"
	else if (is.cqf (object) )
		TRUE
	else
		FALSE
}

is.dks = function (object) inherits (object, "dks")
is.cks = function (object, include.xmix=TRUE) .inc (object, include.xmix, "cks", "xmix")
is.cat = function (object, include.gmix=TRUE) .inc (object, include.gmix, "cat", "gmix")
is.el = function (object) inherits (object, "el")

is.gmix = function (object) inherits (object, "gmix")
is.xmix = function (object) inherits (object, "xmix")

ph.namesf.phmodel = function (sf, ..., all=FALSE)
{	if (is.cat (sf) )
		vars = attr (sf, "gname")
	else
		vars = attr (sf, "xname")

	if (is.gmix (sf) || is.xmix (sf) )
	{	if (all)
			stop ("all needs to be false, for gmix/xmix models")
	}
	else if ( (is.dpdc (sf) || is.cpdc (sf) ) && ! all)
	{	m = attr (sf, "m")
		if (is.pduv (sf) )
			vars = vars [m]
		else
		{	ncon = attr (sf, "ncon")
			vars = vars [(ncon + 1):m]
		}
	}
	vars
}

.range.phpd = function (sf, infv=FALSE, freq=FALSE, n)
{	if (infv)
	{	if (is.cpd (sf) && freq)
			stop ("freq needs to be false, for continuous models")
		if (is.pmf (sf) || is.pdf (sf) )
			c (0, ph.mode (sf, TRUE, freq=freq, n=n) )
		else
		{	if (freq)
			{	if (missing (n) )
				{	if (is.cat (sf) )
						n = sf %$% ".gsum"
					else
						n = sf %$% ".xsum"
				}
				c (0, n)
			}
			else
				c (0, 1)
		}
	}
	else
	{	if (is.cat (sf) )
			attr (sf, "glim")
		else if (is.cks (sf) || is.xmix (sf) )
		{	xlim = attr (sf, "data")$xlim
			if (is.cpdc (sf) )
			{	ncon = attr (sf, "ncon")
				#currently some xmix models based on cksuv, others based on cksc
				if (! is.null (ncon) )
				{	m = attr (sf, "m")
					xlim = xlim [(ncon + 1):m,, drop = is.cpduv (sf)]
				}
			}
			xlim
		}
		else
		{	#kernels
			attr (sf, "xlim")
		}
	}
}

.phminmax = function (sf, infv, which, freq=FALSE, n)
{	k = .range.phpd (sf, infv, freq=freq, n=n)
	if (is.pduv (sf) || infv)
		k [which]
	else
		k [,which]
}

range.dpd = function (sf, infv=FALSE, ..., freq=FALSE, n)
	.range.phpd (sf, infv, freq, n)
range.cpd = function (sf, infv=FALSE, ...)
	.range.phpd (sf, infv)

min.dpd = function (sf, infv=FALSE, ..., freq=FALSE, n)
	.phminmax (sf, infv, 1, freq, n)
max.dpd = function (sf, infv=FALSE, ..., freq=FALSE, n)
	.phminmax (sf, infv, 2, freq, n)
min.cpd = function (sf, infv=FALSE, ...)
	.phminmax (sf, infv, 1)
max.cpd = function (sf, infv=FALSE, ...)
	.phminmax (sf, infv, 2)

seq.dpduv = function (sf, infv=FALSE, ..., midpoints=TRUE, freq=FALSE, n)
{	if (infv)
	{	if (is.pmf (sf) )
		{	if (freq == TRUE)
			{	x = seq (sf)
				sf (x, freq=TRUE, n=n)
			}
			else
				sf %$% ".probs"
		}
		else if (is.dcdf (sf) )
			sf %$% ".PROBS"
		else
		{	P = sf %$% ".PROBS"
			P = sort (unique (P) )
			P = c (0, P)
			if (midpoints)
				P = .midpoints (P)
			P
		}
	}
	else
	{	xlim = range (sf)
		xlim [1]:xlim [2]
	}
}

seq.cpduv = function (sf, infv=FALSE, ..., n=200)
{	if (infv)
	{	if (is.cqf (sf) )
			seq (0, 1, length.out=n)
		else
		{	xlim = range (sf)
			sf (seq (xlim [1], xlim [2], length.out=n) )
		}
	}
	else
	{	xlim = range (sf, infv)
		seq (xlim [1], xlim [2], length.out=n)
	}
}

ph.printf.phmodel = function (sf, ...)
	.object.summary (sf)
ph.printf.dset = function (vf, ...)
	print (paste (class (vf)[1], "object") )

.plot.kernel = function (k, continuous, cdf, ..., main = k %$% "name")
{	if (cdf) sf = k@F
	else sf = k@f
	attr (sf, "variable.name") = "x"
	if (continuous)
	{	attr (sf, "xlim") = c (-1, 1)
		plot_cpd (sf, ..., main=main)
	}
	else
	{	attr (sf, "xlim") = k@xlim
		attr (sf, "freq") = FALSE
		plot_dpd (sf, ..., main=main)
	}
}

ph.plotf.DKernel = function (dk, ..., cdf=FALSE) .plot.kernel (dk, FALSE, cdf, ...)
ph.plotf.CKernel = function (ck, ..., cdf=FALSE) .plot.kernel (ck, TRUE, cdf, ...)

ph.plotf.dksuv = function (sf, data=FALSE, ...) plot_dpd (sf, data, ...)
ph.plotf.catuv = function (sf, ...) plot_dpd (sf, FALSE, ...)
ph.plotf.catc = function (sf, ...) plot_dpd (sf, FALSE, ...)
ph.plotf.gmix = function (sf, ...) plot_dpd (sf, FALSE, ...)

ph.plotf.cksuv = function (sf, data=FALSE, ...)
{	if (is.cqf (sf) || is.cpdc (sf) )
		data = FALSE
	plot_cpd (sf, data, ...)
}

ph.plotf.cksmv = function (sf, in3d=FALSE, data=FALSE,...)
{	if (is.cpdmvc (sf) )
		M = sf %$% "M"
	else
		M = sf %$% "m"
	if (M == 2) plot_cpd_bv (sf, in3d, data, ...)
	else if (M == 3) plot_cpd_tv (sf, in3d, ...)
	else stop ("can only plot mv PDFs/CDFs with 2/3 RVs")
}

ph.plotf.cksc = function (sf, ...) ph.plotf.cksuv (sf, FALSE, ...)
ph.plotf.cksmvc = function (sf, in3d=FALSE, data=FALSE, ...) ph.plotf.cksmv (sf, in3d, data, ...)

ph.plotf.eluv = function (sf, data=FALSE, ...)
{	plot_cpd (sf, FALSE, ...)
	if (data)
	{	x = sf %$% "spline.function" %$% "cx"
		y = sf %$% "spline.function" %$% "cy"
		if (is.cqf (sf) )
			points (x, y, pch=16)
		else
			points (x, y, pch=16)
	}
}

ph.plotf.xmix = function (sf, ...) ph.plotf.cksuv (sf, FALSE, ...)

ph.linesf.cpduv = function (sf, ..., xlim, n=200)
{	if (is.cqf (sf) )
	{	x = seq (0, 1, length.out=n)
		y = sf (x)
	}
	else
	{	if (missing (xlim) )
			xlim = range (sf)
		x = seq (xlim [1], xlim [2], length.out=n)
		y = sf (x)
	}
	lines (x, y, ...)
}

mf.dfh = function (x, ..., freq=FALSE, n) 0
mf.dFh = function (x, ..., freq=FALSE, n) 0
mf.dFht = function (p) 0
mf.gfh = function (g, ..., freq=FALSE, n) 0
mf.gFh = function (g, ..., freq=FALSE, n) 0
mf.gFht = function (p, ..., level.names=FALSE) 0
mf.cfh = function (x) 0
mf.cFh = function (x) 0
mf.cFht = function (p) 0
mf.cfh.mv = function (x) 0
mf.cFh.mv = function (x) 0
mf.chFht = function (p) 0

.arg.error = function (...)
{	expr = format ( (sys.call (-1)) )
	n = length (list (...) )
	if (n > 0)
	{	cat ("call with unsupported args:\n")
		print (expr)
		cat ("check for incorrect argument names\n")
		cat ("check for unnamed non-leading arguments\n")
		cat ("e.g. in pdfuv.cks (x, bw), bw\n")
		warning ("unsupported args, check arg names")
	}
}
