#probhat: Multivariate Generalized Kernel Smoothing and Related Statistical Methods
#Copyright (C), Abby Spurdle, 2018 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

prep.ph.data = function (..., eval=TRUE, echo=FALSE)
{	lines = character (31)
	lines [1] = 'data (Traffic, package="MASS")'
	lines [2] = 'traffic.table <- table (Traffic$y [Traffic$limit=="yes"])'
	lines [3] = 'traffic.bins <- as.integer (names (traffic.table) )'
	lines [4] = 'traffic.bins <- cbind (naccidents=traffic.bins)'
	lines [5] = 'traffic.freq <- as.vector (traffic.table)'
	lines [6] = ''
	lines [7] = 'trees2 <- as.matrix (datasets::trees)[,c (2, 1, 3)]'
	lines [8] = 'colnames (trees2) <- tolower (colnames (trees2) )'
	lines [9] = '#Height (-> m)'
	lines [10] = 'trees2 [,"height"] <- 0.3048 * trees2 [,"height"]'
	lines [11] = '#Girth (-> cm)'
	lines [12] = 'trees2 [,"girth"] <- 2.54 * trees2 [,"girth"]'
	lines [13] = '#Volume (-> m ^ 3)'
	lines [14] = 'trees2 [,"volume"] <- 0.0283168 * trees2 [,"volume"]'
	lines [15] = ''
	lines [16] = 'height <- trees2 [,"height", drop=FALSE]'
	lines [17] = ''
	lines [18] = 'quakes2 <- as.matrix (datasets::quakes)[,c (2, 1, 3:4)]'
	lines [19] = ''
	lines [20] = 'crimes <- cbind (state.x77 [,1, drop=FALSE] * 1e3, as.matrix (USArrests [,-3]) / 1e5)'
	lines [21] = ''
	lines [22] = 'crime.type <- as.factor (rep (colnames (crimes)[-1], each=50) )'
	lines [23] = 'crime.type <- list (crime.type=crime.type)'
	lines [24] = 'n.arrests <- as.vector (round (crimes [,1] * crimes [,-1]) )'
	lines [25] = ''
	lines [26] = 'species <- list (species=iris$Species)'
	lines [27] = 'sepal.length <- cbind (sepal.length = iris$Sepal.Length)'
	lines [28] = 'sepal.width <- cbind (sepal.width = iris$Sepal.Width)'
	lines [29] = ''
	lines [30] = 'data (unemployment, package="fclust")'
	lines [31] = 'unemployment <- as.matrix (unemployment)[,-2]'
	if (eval)
		eval (parse (text=lines), .GlobalEnv)
	if (echo)
	{	for (line in lines)
			cat (line, "\n", sep="")
	}
}
