#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.n.unique = function (x) length (unique (x) )
.any.duplicates = function (x) .n.unique (x) != length (x)

.seq.dks = function (lim) lim [1]:lim [2]
.seq.cks = function (lim, n) seq (lim [1], lim [2], length.out=n)

range.pmfuv = function (f, ...) f %$% "xlim"
range.pdfuv = function (f, ...) f %$% "xlim"
range.cdfuv = function (F, ...) F %$% "xlim"
range.qfuv = function (F.inv, ...) F.inv %$% "xlim"

seq.pmfuv = function (f, ...) .seq.dks (range (f) )
seq.pdfuv = function (f, n, ...) .seq.cks (range (f), n)

seq.cdfuv = function (F, n, ...)
{	if (is.dpd (F) )
		.seq.dks (range (F) )
	else
		.seq.cks (range (F), n)
}

seq.qfuv = function (F.inv, n, ...)
	seq.cdfuv (F.inv, n)
