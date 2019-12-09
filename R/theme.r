#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

use.ph.theme = function (theme="blue",
	area.color, area.color.2="grey75", palette="Earth")
{	pho = list ()
	if (theme == "blue")
	{	use.theme ("blue")
		bso = getOption ("barsurf")
		bso$plot2d.contour.colv.1 = c (157.5, 30, 80)
		bso$plot2d.contour.colv.2 = c (282.5, 30, 80)
		options (barsurf=bso)
		if (missing (area.color) )
			area.color = rgb (0.3, 0.6, 1, 0.55)
	}
	else if (theme == "green")
	{	use.theme ("green")
		bso = getOption ("barsurf")
		bso$plot2d.contour.colv.1 = c (200, 20, 80)
		options (barsurf=bso)
		if (missing (area.color) )
			area.color = rgb (0, 0.55, 0, 0.55)
	}
	else
		stop ("unsupported theme")
	pho$area.color = area.color
	pho$area.color.2 = area.color.2
	pho$palette = palette
	options (probhat=pho)
}
