### R code from vignette source 'probhat.Rnw'

###################################################
### code chunk number 1: probhat.Rnw:35-39
###################################################
options(continue=" ", width=80)
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 2.6, 1.6), cex=0.7, cex.main=1)))
set.seed (1)


###################################################
### code chunk number 2: probhat.Rnw:152-155
###################################################
library (probhat)       #required
library (fclust)        #optional, used in appendix
library (scatterplot3d) #optional


###################################################
### code chunk number 3: probhat.Rnw:162-163
###################################################
set.ph.options (rendering.style="pdf", theme="green")


###################################################
### code chunk number 4: probhat.Rnw:168-169
###################################################
prep.ph.data ()


###################################################
### code chunk number 5: probhat.Rnw:183-186
###################################################
dfh <- pmfuv.dks (traffic.bins, traffic.freq)
dFh <- cdfuv.dks (traffic.bins, traffic.freq)
dFht <- qfuv.dks (traffic.bins, traffic.freq)


###################################################
### code chunk number 6: probhat.Rnw:190-191
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dfh, TRUE, main="Probability Mass Function")


###################################################
### code chunk number 7: probhat.Rnw:193-194
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dfh, TRUE, freq=TRUE, main="same as above\nbut with frequencies")


###################################################
### code chunk number 8: probhat.Rnw:196-197
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dFh, TRUE, main="Cumulative Distribution Function")


###################################################
### code chunk number 9: probhat.Rnw:199-200
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dFht, main="Quantile Function")


###################################################
### code chunk number 10: probhat.Rnw:209-212
###################################################
dfh (10)
dfh (10, freq=TRUE)
dFht (c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 11: probhat.Rnw:230-233
###################################################
cfh <- pdfuv.cks (height)
cFh <- cdfuv.cks (height)
cFht <- qfuv.cks (height)


###################################################
### code chunk number 12: probhat.Rnw:237-238
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh, TRUE, main="Probability Density Function")


###################################################
### code chunk number 13: probhat.Rnw:240-241
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFh, TRUE, main="Cumulative Distribution Function")


###################################################
### code chunk number 14: probhat.Rnw:243-244
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFht, main="Quantile Function")


###################################################
### code chunk number 15: probhat.Rnw:251-253
###################################################
cfh (22)
cFht (c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 16: probhat.Rnw:264-267
###################################################
#note, depth is third variable
#(will fit a semi-bounded model)
head (quakes2, 2)


###################################################
### code chunk number 17: probhat.Rnw:270-273
###################################################
cfh2 <- pdfmv.cks (quakes2 [,1:2], smoothness = c (0.35, 1) )
cfh3 <- pdfmv.cks (quakes2 [,1:3], smoothness = c (0.35, 1, 1),
    a = c (-Inf, -Inf, 0) )


###################################################
### code chunk number 18: probhat.Rnw:277-279
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2,, TRUE,
    main="Bivariate PDF, 2D")


###################################################
### code chunk number 19: probhat.Rnw:281-282
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2, TRUE, main="Bivariate PDF, 3D")


###################################################
### code chunk number 20: probhat.Rnw:284-286
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh3, main="Trivariate PDF",
    nslides=8, zlim = c (680, 40) )


###################################################
### code chunk number 21: probhat.Rnw:294-296
###################################################
cfh2 (c (180, -20) )
cfh3 (c (180, -20, 300) )


###################################################
### code chunk number 22: probhat.Rnw:312-313
###################################################
conditions <- c (long=180, lat=-20)


###################################################
### code chunk number 23: probhat.Rnw:316-321
###################################################
depth.fhc <- pdfc.cks (quakes2 [,-4], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE,
    a = c (-Inf, -Inf, 0) )
mag.fhc <- pdfc.cks (quakes2 [,-3], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE)


###################################################
### code chunk number 24: probhat.Rnw:325-326
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.fhc, main="conditional distribution of depth\n(long=180, lat=-20)")


###################################################
### code chunk number 25: probhat.Rnw:328-329
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (mag.fhc, main="conditional distribution of mag\n(long=180, lat=-20)")


###################################################
### code chunk number 26: probhat.Rnw:342-348
###################################################
depth.mag.fhc <- pdfmvc.cks (quakes2, smoothness = c (0.35, 1, 1, 1),
    conditions = c (long=180, lat=-20), preserve.range=TRUE,
    a = c (-Inf, -Inf, 0, -Inf) )
lat.long.fhc <- pdfmvc.cks (quakes2 [,-4], smoothness = c (0.35, 1, 1),
    conditions = c (depth=168), preserve.range=TRUE,
    a = c (-Inf, -Inf, 0) )


###################################################
### code chunk number 27: probhat.Rnw:352-354
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.mag.fhc,
    main="conditional distribution of depth and mag\n(long=180, lat=-20)")


###################################################
### code chunk number 28: probhat.Rnw:356-358
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (lat.long.fhc,
    main="conditional distribution of lat and long\n(depth=168)")


###################################################
### code chunk number 29: probhat.Rnw:372-373
###################################################
gfh <- pmfuv.cat (crime.type, n.arrests)


###################################################
### code chunk number 30: probhat.Rnw:377-378
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, main="Probability Mass Function")


###################################################
### code chunk number 31: probhat.Rnw:380-381
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, freq=TRUE, main="same as above\nbut with frequencies")


###################################################
### code chunk number 32: probhat.Rnw:388-389
###################################################
levels (crime.type$crime.type)


###################################################
### code chunk number 33: probhat.Rnw:392-394
###################################################
gfh (1)
gfh ("Assault")


###################################################
### code chunk number 34: probhat.Rnw:397-398
###################################################
gfh ("Assault", freq=TRUE)


###################################################
### code chunk number 35: probhat.Rnw:414-416
###################################################
eFh <- cdf.el (height)
eFht <- qf.el (height)


###################################################
### code chunk number 36: probhat.Rnw:420-421
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFh, main="Cumulative Distribution Function")


###################################################
### code chunk number 37: probhat.Rnw:423-424
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFht, main="Quantile Function")


###################################################
### code chunk number 38: probhat.Rnw:451-455
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1a.gmix <- ph4.pmfc.gmix (species, sepal.length,
    conditions = c (sepal.length=5.5) )
plot (fh1a.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=5.5)")


###################################################
### code chunk number 39: probhat.Rnw:457-461
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1b.gmix <- ph4.pmfc.gmix (species, sepal.length,
    conditions = c (sepal.length=6.5) )
plot (fh1b.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=6.5)")


###################################################
### code chunk number 40: probhat.Rnw:469-477
###################################################
getOption("SweaveHooks")[["fig"]]()
fh2.gmix <- ph4.pmfc.gmix (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )
plot (fh2.gmix,
    main = paste (
        "Conditional Distribution of Iris Species",
        "(sepal.length=6, sepal.width=3)",
        sep="\n")
    )


###################################################
### code chunk number 41: probhat.Rnw:486-488
###################################################
ph.mode (fh2.gmix)
ph.mode (fh2.gmix, level.names=TRUE)


###################################################
### code chunk number 42: probhat.Rnw:508-509
###################################################
gset <- ph4.pdfuv.gset.cks (species, sepal.length)


###################################################
### code chunk number 43: probhat.Rnw:513-514
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gset, main="Density Estimates of Sepal Length\n(grouped by species)")


###################################################
### code chunk number 44: probhat.Rnw:520-521
###################################################
mset = ph4.qfuv.mset.el (trees2)


###################################################
### code chunk number 45: probhat.Rnw:525-526
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (mset, nr=2, nc=2)


###################################################
### code chunk number 46: probhat.Rnw:545-548
###################################################
#multivariate cdf
#note pwith function can use pdf or cdf, but only nonconditional
cFh3 <- cdfmv.cks (trees2)


###################################################
### code chunk number 47: probhat.Rnw:551-557
###################################################
xlim <- matrix (c (
    22, 24,    #height in 22 to 24
    28, 38,    #girth  in 28 to 38
    0.55, 1.05 #volume in 0.55 to 1.05
    ),, 2, byrow=TRUE, dimnames = list (colnames (trees2), c ("a", "b") ) )
xlim


###################################################
### code chunk number 48: probhat.Rnw:560-563
###################################################
#multivariate probability
probmv (cFh3, xlim [,1], xlim [,2])
pwith (cFh3, xlim=xlim)


###################################################
### code chunk number 49: probhat.Rnw:574-580
###################################################
gset1a <- ph4.pdfuv.gset.cks (species, sepal.length)
gset1b <- ph4.pdfuv.gset.cks (species, sepal.width)
gset2 <- ph4.pdfmv.gset.cks (species, cbind (sepal.length, sepal.width) )
dists1a <- pdist (gset1a)
dists1b <- pdist (gset1b)
dists2 <- pdist (gset2)


###################################################
### code chunk number 50: probhat.Rnw:583-589
###################################################
names <- names (gset2)
cols <- c ("red", "blue", "darkgreen")
main <- paste0 (substring (names, 1, 4), " (", cols, ")", collapse=", ")
I1 <- (names [1] == species [[1]])
I2 <- (names [2] == species [[1]])
I3 <- (names [3] == species [[1]])


###################################################
### code chunk number 51: probhat.Rnw:593-604
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gset2 [[1]], contour.color = cols [1],
    main=main,
    xlim = range (sepal.length), ylim = range (sepal.width),
    heatmap=FALSE)
plot (gset2 [[2]], contour.color = cols [2],
    add=TRUE, heatmap=FALSE)
plot (gset2 [[3]], contour.color = cols [3],
    add=TRUE, heatmap=FALSE)
points (sepal.length [I1], sepal.width [I1], pch=16, col = cols [1])
points (sepal.length [I2], sepal.width [I2], pch=16, col = cols [2])
points (sepal.length [I3], sepal.width [I3], pch=16, col = cols [3])


###################################################
### code chunk number 52: probhat.Rnw:612-620
###################################################
#distance matrix
#(both variables)
dists2
#pairs, ranked
ph4.rdist (dists2)  #both
ph4.rdist (dists1a) #length
ph4.rdist (dists1b) #width



###################################################
### code chunk number 53: probhat.Rnw:644-646
###################################################
chFht <- chqf.cks (trees2)
synthetic.data <- rng (chFht, 31)


###################################################
### code chunk number 54: probhat.Rnw:649-657
###################################################
#convenience function
plot.trees.data <- function (x, main)
{   height <- x [,"height"]
    girth <- x [,"girth"]
    volume <- x [,"volume"]
    scatterplot3d (height, girth, volume,
        main=main, type="h", angle=112.5, pch=16)
}


###################################################
### code chunk number 55: probhat.Rnw:661-663
###################################################
getOption("SweaveHooks")[["fig"]]()
#original data
plot.trees.data (trees2, "original data")


###################################################
### code chunk number 56: probhat.Rnw:666-668
###################################################
getOption("SweaveHooks")[["fig"]]()
#synthetic data
plot.trees.data (synthetic.data, "synthetic data")


###################################################
### code chunk number 57: probhat.Rnw:705-709
###################################################
ph.mean (cFh)
ph.var (cFh)
ph.skewness (cFh)
ph.kurtosis (cFh)


###################################################
### code chunk number 58: probhat.Rnw:721-723
###################################################
quartiles (cFht)
deciles (cFht)


###################################################
### code chunk number 59: probhat.Rnw:726-728
###################################################
ntiles (8, cFht)
ntiles (8, cFht, prob=TRUE)


###################################################
### code chunk number 60: probhat.Rnw:740-742
###################################################
ph.median (cFht)
ph.quantile (cFht, c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 61: probhat.Rnw:745-747
###################################################
#inter-quartile range
iqr (cFht)


###################################################
### code chunk number 62: probhat.Rnw:750-752
###################################################
#inter-quantile ranges
iqr (cFht, 2/3)


###################################################
### code chunk number 63: probhat.Rnw:760-762
###################################################
cFht (0.5)
ph.median (cFht)


###################################################
### code chunk number 64: probhat.Rnw:766-767
###################################################
ph.median (height)


###################################################
### code chunk number 65: probhat.Rnw:779-781
###################################################
ph.mode(cfh)
ph.mode(cfh, TRUE)


###################################################
### code chunk number 66: probhat.Rnw:788-797
###################################################
height.summary <- c (
    mean = ph.mean (cFh),         #CDF
    sd = ph.sd (cFh),             #CDF
    variance = ph.var (cFh),      #CDF
    skewness = ph.skewness (cFh), #CDF
    kurtosis = ph.kurtosis (cFh), #CDF
                                  #
    median = ph.median (cFht),    #QF
    mode = ph.mode (cfh) )        #PDF


###################################################
### code chunk number 67: probhat.Rnw:800-805
###################################################
strs <- c (c ("mean", "median", "mode") )
x <- height.summary [strs]
y <- c (0.06, 0.1, 0.14)
colors <- c ("black", "blue", "darkgreen")
adjv <- c (1.25, 0.5, -0.25)


###################################################
### code chunk number 68: probhat.Rnw:809-813
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh)
abline (v=x, col=colors)
for (i in 1:3)
    text (x [i], y [i], strs [i], adj = adjv [i], col = colors [i])


###################################################
### code chunk number 69: probhat.Rnw:817-818
###################################################
height.summary


###################################################
### code chunk number 70: probhat.Rnw:1147-1148
###################################################
membership <- FKM.gk (unemployment, k=3, seed=2)$U


###################################################
### code chunk number 71: probhat.Rnw:1153-1155
###################################################
w <- membership [,1]
w <- w / sum (w)


###################################################
### code chunk number 72: probhat.Rnw:1160-1161
###################################################
wfh.1 <- pdfmv.cks (unemployment, w=w)


###################################################
### code chunk number 73: probhat.Rnw:1165-1166
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (wfh.1,, TRUE)


###################################################
### code chunk number 74: probhat.Rnw:1169-1171
###################################################
getOption("SweaveHooks")[["fig"]]()
k = 1 - w / max (w)
plot (unemployment, pch=16, col=rgb (k, k, k) )


###################################################
### code chunk number 75: probhat.Rnw:1176-1180
###################################################
w <- membership [,2]
wfh.2 = pdfmv.cks (unemployment, w = w / sum (w) )
w <- membership [,3]
wfh.3 = pdfmv.cks (unemployment, w = w / sum (w) )


###################################################
### code chunk number 76: probhat.Rnw:1184-1187 (eval = FALSE)
###################################################
## plot (wfh.1, main="cluster 1")
## plot (wfh.2, main="cluster 2")
## plot (wfh.3, main="cluster 3")


###################################################
### code chunk number 77: probhat.Rnw:1191-1197
###################################################
getOption("SweaveHooks")[["fig"]]()
p0 <- par (mfrow = c (2, 2) )
plot (wfh.1, main="cluster 1")
plot.new ()
plot (wfh.2, main="cluster 2")
plot (wfh.3, main="cluster 3")
par (p0)


###################################################
### code chunk number 78: probhat.Rnw:1205-1206
###################################################
prep.ph.data (eval=FALSE, echo=TRUE)


###################################################
### code chunk number 79: probhat.Rnw:1214-1218
###################################################
ht <- function (object)
{   print (head (object, 2) )
    print (tail (object, 2) )
}


###################################################
### code chunk number 80: probhat.Rnw:1221-1222
###################################################
ht (cbind (traffic.bins, freq=traffic.freq) )


###################################################
### code chunk number 81: probhat.Rnw:1225-1231
###################################################
ht (trees2)
ht (quakes2)
ht (crimes)
ht (data.frame (crime.type=crime.type$crime.type, n.arrests) )
ht (data.frame (species, sepal.length, sepal.width) )
ht (unemployment)


