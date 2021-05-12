#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

ph.namesf = function (...) UseMethod ("ph.namesf")
ph.printf = function (...) UseMethod ("ph.printf")
ph.plotf = function (...) UseMethod ("ph.plotf")
ph.linesf = function (...) UseMethod ("ph.linesf")

names.phob = function (x, ...) ph.namesf (x, ...)
print.phob = function (x, ...) ph.printf (x, ...)
plot.phob = function (x, ...) ph.plotf (x, ...)
lines.phob = function (x, ...) ph.linesf (x, ...)
