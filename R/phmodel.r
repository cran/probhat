#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

print.phmodel = function (x, ...)
	object.summary (x, ...)

is.dpd = function (f) inherits (f, "dpd")
is.cpd = function (f) inherits (f, "cpd")
is.pmfuv = function (f) inherits (f, "pmfuv")
is.pdfuv = function (f) inherits (f, "pdfuv")
is.cdfuv = function (f) inherits (f, "cdfuv")
is.qfuv = function (f) inherits (f, "qfuv")
is.chqf = function (f) inherits (f, "chqf")
