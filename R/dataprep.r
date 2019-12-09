#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

data.prep = function (eval=TRUE, echo=FALSE)
{	lines = character (30)
	lines [1] = 'data (Traffic, package="MASS")'
	lines [2] = 'traffic.table = table (Traffic$y [Traffic$limit=="yes"])'
	lines [3] = 'traffic.x = as.integer (names (traffic.table) )'
	lines [4] = 'traffic.x = COL (traffic.x, "naccidents")'
	lines [5] = 'traffic.h = as.vector (traffic.table)'
	lines [6] = ''
	lines [7] = 'region = COL (datasets::state.region, "region")'
	lines [8] = ''
	lines [9] = 'iris = datasets::iris'
	lines [10] = 'iris.Species = as.character (iris$Species)'
	lines [11] = 'iris.Species = COL (iris.Species, "Species")'
	lines [12] = 'iris.Sepal.Length = iris$Sepal.Length'
	lines [13] = 'iris.Sepal.Length = COL (iris.Sepal.Length, "Sepal Length")'
	lines [14] = ''
	lines [15] = 'quakes = datasets::quakes'
	lines [16] = 'quakes = as.matrix (quakes)[,-5]'
	lines [17] = ''
	lines [18] = 'trees = datasets::trees'
	lines [19] = 'trees = as.matrix (trees)'
	lines [20] = '#Girth (-> cm)'
	lines [21] = 'trees [,"Girth"] = 2.54 * trees [,"Girth"]'
	lines [22] = '#Height (-> m)'
	lines [23] = 'trees [,"Height"] = 0.3048 * trees [,"Height"]'
	lines [24] = '#Volume (-> m ^ 3)'
	lines [25] = 'trees [,"Volume"] = 0.0283168 * trees [,"Volume"]'
	lines [26] = ''
	lines [27] = 'Height = COL.of (trees, "Height")'
	lines [28] = ''
	lines [29] = 'data (unemployment, package="fclust")'
	lines [30] = 'unemployment = as.matrix (unemployment)[,-2]'
	if (eval)
	{	for (line in lines)
			eval (parse (text=line), .GlobalEnv)
	}
	if (echo)
	{	for (line in lines)
			cat (line, "\n", sep="")
	}
}
