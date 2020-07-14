#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.onLoad = function (...)
	set.ph.options ()

set.ph.options = function (...,
	rendering.style="r", theme="blue",
	main.line.width=1, main.line.color="#000000", main.fill.color,
	main.fill.color.2="#B0B0B0", semi.fill.color="#00000030",
	palette="earth")
{	set.bs.options (rendering.style=rendering.style, theme=theme)
	if (missing (main.fill.color) )
	{	if (theme == "blue")
			main.fill.color = rgb (0.3, 0.6, 1, 0.55)
		else if (theme == "green")
			main.fill.color = rgb (0, 0.55, 0, 0.55)
		else
			main.fill.color = main.fill.color.2
	}
	pho = list ()
	pho$main.line.width = main.line.width
	pho$main.line.color = main.line.color
	pho$main.fill.color = main.fill.color
	pho$main.fill.color.2 = main.fill.color.2
	pho$semi.fill.color = semi.fill.color
	pho$palette = palette
	options (probhat=pho)
}

set.ph.theme = function (theme="blue")
{	set.bs.options (theme=theme)
	if (theme == "blue")
		fill.color = rgb (0.3, 0.6, 1, 0.55)
	else if (theme == "green")
		fill.color = rgb (0, 0.55, 0, 0.55)
	else
		stop ("set.ph.theme only supports blue and green")
	pho = getOption ("probhat")
	pho$fill.color = fill.color
	options (probhat=pho)
}
