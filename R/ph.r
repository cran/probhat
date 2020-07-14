#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.KV.pmf = c ("pmfuv", "dpduv", "pmf", "dpd", "phob")
.KV.dcdf = c ("dcdfuv", "dpduv", "dcdf", "dpd", "phob")
.KV.pdf = c ("pcfuv", "cpduv", "pdf", "cpd", "phob")
.KV.ccdf = c ("ccdfuv", "cpduv", "ccdf", "cpd", "phob")

.CV = c ("phmodel", "phob")

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
.CV.pmfc.cat = c ("pmfc.cat", "catc", "catuv", "cat", "pmfc", "pmfuv", "dpdc", "dpduv", "pmf", "dpd", .CV)
.CV.cdfc.cat = c ("cdfc.cat", "catc", "catuv", "cat", "dcdfc", "dcdfuv", "dpdc", "dpduv", "dcdf", "dpd", .CV)
.CV.qfc.cat = c ("qfc.cat", "catc", "catuv", "cat", "dqfc", "dqfuv", "dpdc", "dpduv", "dqf", "dpd", .CV)

.CV.cdfuv.el = c ("cdfuv.el", "eluv", "el", "ccdfuv", "cpduv", "ccdf", "cpd", .CV)
.CV.qfuv.el = c ("qfuv.el", "eluv", "el", "cqfuv", "cpduv", "cqf", "cpd", .CV)

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
dfreq = function (...) pmfuv.dks (..., freq=TRUE)
gfreq = function (...) pmfuv.cat (..., freq=TRUE)

is.phob = function (xf) inherits (xf, "phob")
is.phmodel = function (xf) inherits (xf, "phmodel")

is.dpd = function (xf) inherits (xf, "dpd")
is.cpd = function (xf) inherits (xf, "cpd")

is.pmf = function (xf) inherits (xf, "pmf")
is.dcdf = function (xf) inherits (xf, "dcdf")
is.dqf = function (xf) inherits (xf, "dqf")
is.pdf = function (xf) inherits (xf, "pdf")
is.ccdf = function (xf) inherits (xf, "ccdf")
is.cqf = function (xf) inherits (xf, "cqf")
is.cchqf = function (xf) inherits (xf, "cchqf")

is.dpduv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "dpduv", "dpdc")
is.dpdc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "dpdc", "dpdmvc")
is.cpduv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "cpduv", "cpdc")
is.cpdmv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "cpdmv", "cpdmvc")
is.cpdc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "cpdc", "cpdmvc")
is.cpdmvc = function (xf) inherits (xf, "cpdmvc")

is.pmfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "pmfuv", "pmfc")
is.pmfc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "pmfc", "pmfmvc")
is.dcdfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "dcdfuv", "dcdfc")
is.dcdfc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "dcdfc", "dcdfmvc")
is.dqfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "dqfuv", "dqfc")
is.dqfc = function (xf) inherits (xf, "dqfc")
is.pdfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "pdfuv", "pdfc")
is.pdfmv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "pdfmv", "pdfmvc")
is.pdfc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "pdfc", "pdfmvc")
is.pdfmvc = function (xf) inherits (xf, "pdfmvc")
is.ccdfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "ccdfuv", "ccdfc")
is.ccdfmv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "ccdfmv", "ccdfmvc")
is.ccdfc = function (xf, include.multivariate=TRUE) .inc (xf, include.multivariate, "ccdfc", "ccdfmvc")
is.ccdfmvc = function (xf) inherits (xf, "ccdfmvc")
is.cqfuv = function (xf, include.conditional=TRUE) .inc (xf, include.conditional, "cqfuv", "cqfc")
is.cqfc = function (xf) inherits (xf, "cqfc")

is.phspline = function (xf)
{	if (is.pdfuv (xf, TRUE) || is.ccdfuv (xf, TRUE) )
		xf %$% "is.spline"
	else if (is.cqf (xf) )
		TRUE
	else
		FALSE
}

is.dks = function (xf) inherits (xf, "dks")
is.cks = function (xf) inherits (xf, "cks")
is.cat = function (xf) inherits (xf, "cat")
is.el = function (xf) inherits (xf, "el")

.frange = function (sf, infv, freq)
{	if (infv)
	{	if (is.pmf (sf) || is.pdf (sf) )
			c (0, ph.mode (sf, TRUE, freq=freq) )
		else
			c (0, 1)
	}
	else
		sf %$% "xlim"
}

min.dpduv = function (sf, infv=FALSE, ...) range (sf, infv, ...) [1]
max.dpduv = function (sf, infv=FALSE, ...) range (sf, infv, ...) [2]
min.cpduv = function (sf, infv=FALSE, ...) range (sf, infv, ...) [1]
max.cpduv = function (sf, infv=FALSE, ...) range (sf, infv, ...) [2]

range.dpduv = function (sf, infv=FALSE, ..., freq) .frange (sf, infv, freq)
range.cpduv = function (sf, infv=FALSE, ...) .frange (sf, infv)

seq.dpduv = function (sf, infv=FALSE, ..., midpoints=TRUE, freq)
{	if (infv)
	{	if (is.pmf (sf) )
		{	if (missing (freq) )
				freq = sf %$% "freq"
			if (freq == TRUE)
			{	x = seq (sf)
				sf (x, freq=TRUE)
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

print.phmodel = function (x, ...)
	object.summary (x, ...)

.plot.kernel = function (kc, continuous, cdf, ..., main = kc %$% "name")
{	if (cdf) fs = kc [[2]]
	else fs = kc [[1]]
	if (continuous) plot_cpd (fs, ..., main=main)
	else plot_dpd (fs, ..., main=main)
}

plot.dkernel = function (x, ..., cdf=FALSE) .plot.kernel (x, FALSE, cdf, ...)
plot.ckernel = function (x, ..., cdf=FALSE) .plot.kernel (x, TRUE, cdf, ...)

.plot.dksuv = function (sf, data, ..., combine, space, freq)
{	if (is.dqf (sf) )
		data = FALSE
	if (missing (combine) )
		combine = is.dks (sf)
	if (missing (space) )
	{	if (is.pmf (sf) && is.cat (sf) )
			space = 2
		else
			space = 0
	}
	plot_dpd (sf, data, ..., combine=combine, freq=freq, space=space)
}

plot.dksuv = function (x, data=FALSE, ..., freq)
	.plot.dksuv (x, data, ..., freq=freq)

plot.cksuv = function (x, data=FALSE, ...)
{	sf = x

	if (is.cqf (sf) || is.cpdc (sf) )
		data = FALSE
	plot_cpd (sf, data, ...)
}

plot.cksmv = function (x, in3d=FALSE, data=FALSE,...)
{	sf = x

	if (is.cpdmvc (sf) )
		M = sf %$% "M"
	else
		M = sf %$% "m"
	if (M == 2) plot_cpd_bv (sf, in3d, data, ...)
	else if (M == 3) plot_cpd_tv (sf, ...)
	else stop ("can only plot mv PDFs/CDFs with 2/3 RVs")
}

plot.cksc = function (x, ...) plot.cksuv (x, FALSE, ...)
plot.cksmvc = function (x, in3d=FALSE, data=FALSE, ...) plot.cksmv (x, in3d, data, ...)

plot.catuv = function (x, ..., combine, freq, space)
	.plot.dksuv (x, FALSE, ..., combine=combine, freq=freq, space=space)

plot.eluv = function (x, data=FALSE, ...)
{	sf = x

	plot_cpd (sf, FALSE, ...)
	if (data)
	{	x = sf %$% "spline.function" %$% "cx"
		y = sf %$% "spline.function" %$% "cy"
		if (is.cqf (sf) )
			points (x, y, pch=16)
		else
			points (x, y, pch=16)
	}
}

lines.cpduv = function (x, ...)
{	N = 200
	sf = x

	if (is.cqf (sf) )
	{	x = seq (0, 1, length.out=N)
		y = sf (x)
	}
	else
	{	x = seq (sf, n=N)
		y = sf (x)
	}
	lines (x, y, ...)
}

mf.dfh = function (x, ..., freq) 0
mf.dFh = function (q) 0
mf.dFht = function (p) 0
mf.gfh = function (g, ..., freq) 0
mf.gFh = function (q) 0
mf.gFht = function (p, ..., name=FALSE) 0
mf.cfh = function (x) 0
mf.cFh = function (q) 0
mf.cFht = function (p) 0
mf.cfh.mv = function (x) 0
mf.cFh.mv = function (q) 0
mf.chFht = function (p) 0
