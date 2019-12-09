#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.deflab = function (f, lab)
{	if (missing (lab) )
	{	if (any (class (f) %in% c ("pdfc", "cdfc", "qfc") ) )
			(f %$% "variable.names")[f %$% "m"]
		else
			f %$% "variable.name"
	}
	else
		lab
}

plot.pmfuv.dks = function (x, with.data.bars=FALSE,
	main, xlab, ylab="mass", line.width=1, line.color="black", area.color, ...)
{	f = x

	xlab = .deflab (f, xlab)
	if (with.data.bars)
		.plot.dksuv.2 (f, main, xlab, ylab, line.width, line.color, area.color)
	else
		.plot.dksuv (f, main, xlab, ylab, line.width, line.color, area.color, ...)
}

plot.cdfuv.dks = function (x, with.data.bars=FALSE,
	main, xlab, ylab="cumprob", line.width=1, line.color="black", area.color, ...)
{	F = x

	xlab = .deflab (F, xlab)
	if (with.data.bars)
		.plot.dksuv.2 (F, main, xlab, ylab, line.width, line.color, area.color)
	else
		.plot.dksuv (F, main, xlab, ylab, line.width, line.color, area.color, ...)
}

plot.qfuv.dks = function (x, no.data.bars,
	main, xlab="cumprob", ylab, line.width=1, line.color="black", area.color, ...)
{	F.inv = x

	ylab = .deflab (F.inv, ylab)
	.plot.dksuv.qf (F.inv, main, xlab, ylab, line.width, line.color, area.color, ...)
}

plot.pdfuv.cks = function (x, with.data.points=FALSE,
	main, xlab, ylab="density", line.width=1, line.color="black", area.color, ...,
	add=FALSE)
{	f = x

	xlab = .deflab (f, xlab)
	if (add)
		.add.plot.cksuv (f, line.width, line.color, area.color)
	else if (with.data.points)
		.plot.cksuv.2 (f, main, xlab, ylab, line.width, line.color, area.color)
	else
		.plot.cksuv (f, main, xlab, ylab, line.width, line.color, area.color, ...)
}

plot.cdfuv.cks = function (x, with.data.points=FALSE,
	main, xlab, ylab="cumprob", line.width=1, line.color="black", area.color, ...,
	add=FALSE)
{	F = x

	xlab = .deflab (F, xlab)
	if (add)
		.add.plot.cksuv (F, line.width, line.color, area.color)
	else if (with.data.points)
		.plot.cksuv.2 (F, main, xlab, ylab, line.width, line.color, area.color)
	else
		.plot.cksuv (F, main, xlab, ylab, line.width, line.color, area.color, ...)
}

plot.qfuv.cks = function (x, no.data.points,
	main, xlab="cumprob", ylab, line.width=1, line.color="black", area.color, ...)
{	F.inv = x

	ylab = .deflab (F.inv, ylab)
	.plot.cksuv (F.inv, main, xlab, ylab, line.width, line.color, area.color, ..., is.qf=TRUE)
}

plot.pdfc.cks = function (x, no.data.points,
	main, xlab, ylab="density", line.width=1, line.color="black", area.color, ...,
	add=FALSE)
	plot.pdfuv.cks (x, FALSE, main, xlab, ylab, line.width, line.color, area.color, ..., add=add)

plot.cdfc.cks = function (x, no.data.points,
	main, xlab, ylab="cumprob", line.width=1, line.color="black", area.color, ...,
	add=FALSE)
	plot.cdfuv.cks (x, FALSE, main, xlab, ylab, line.width, line.color, area.color, ..., add=add)

plot.qfc.cks = function (x, no.data.points,
	main, xlab="cumprob", ylab, line.width=1, line.color="black", area.color, ...)
	plot.qfuv.cks (x, FALSE, main, xlab, ylab, line.width, line.color, area.color, ...)

plot.pmfuv.cat = function (x, main, xlab, ylab="mass", line.width=1, line.color="black", area.color, ...)
	.plot.categorical (x, main, .deflab (x, xlab), ylab, line.width, line.color, area.color)

plot.cdfuv.cat = function (x, main, xlab, ylab="cumprob", line.width=1, line.color="black", area.color, ...)
	.plot.categorical (x, main, .deflab (x, xlab), ylab, line.width, line.color, area.color)

plot.qfuv.cat = function (x, main, xlab="cumprob", ylab, line.width=1, line.color="black", area.color, ...)
	.plot.categorical.qf (x, main, xlab, .deflab (x, ylab), line.width, line.color, area.color)

plot.cdf.el = function (x, with.data.points=FALSE,
	main, xlab, ylab="cumprob", line.width=1, line.color="black", area.color, ...)
{	F = x

	xlab = .deflab (F, xlab)
	plot.cdfuv.cks (F, FALSE, main, xlab, ylab, line.width, line.color, area.color)
	if (with.data.points)
	{	x = F %$% "x"
		y = F %$% "spline.function" %$% "cy"
		points (x, y, pch=16)
	}
}

plot.qf.el = function (x, with.data.points=FALSE,
	main, xlab="cumprob", ylab, line.width=1, line.color="black", area.color, ...)
{	F.inv = x

	ylab = .deflab (F.inv, ylab)
	plot.qfuv.cks (F.inv, FALSE, main, xlab, ylab, line.width, line.color, area.color)
	if (with.data.points)
	{	x = F.inv %$% "x"
		y = F.inv %$% "spline.function" %$% "cx"
		points (x, y, pch=16)
	}
}

.ph.lines = function (f, ...)
{	x = seq (f, 200)
	y = f (x)
	lines (x, y, ...)
}

.ph.lines.2 = function (f, ...)
{	x = seq (0, 1, length.out=200)
	y = f (x)
	lines (x, y, ...)
}

lines.pdfuv = function (x, ...) .ph.lines (x, ...)
lines.cdfuv = function (x, ...) .ph.lines (x, ...)
lines.qfuv = function (x, ...) .ph.lines.2 (x, ...)
