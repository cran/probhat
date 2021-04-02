### R code from vignette source 'probhat.Rnw'

###################################################
### code chunk number 1: probhat.Rnw:35-39
###################################################
options(continue=" ", width=80)
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 2.6, 1.6), cex=0.7, cex.main=1)))
set.seed (1)


###################################################
### code chunk number 2: probhat.Rnw:154-157
###################################################
library (probhat)       #required
library (fclust)        #optional, used in appendix
library (scatterplot3d) #optional


###################################################
### code chunk number 3: probhat.Rnw:164-165
###################################################
set.ph.options (rendering.style="pdf", theme="green")


###################################################
### code chunk number 4: probhat.Rnw:170-171
###################################################
prep.ph.data ()


###################################################
### code chunk number 5: probhat.Rnw:185-188
###################################################
dfh <- pmfuv.dks (traffic.bins, traffic.freq)
dFh <- cdfuv.dks (traffic.bins, traffic.freq)
dFht <- qfuv.dks (traffic.bins, traffic.freq)


###################################################
### code chunk number 6: probhat.Rnw:192-193
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dfh, TRUE, main="Probability Mass Function")


###################################################
### code chunk number 7: probhat.Rnw:195-196
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dfh, TRUE, freq=TRUE, main="same as above\nbut with frequencies")


###################################################
### code chunk number 8: probhat.Rnw:198-199
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dFh, TRUE, main="Cumulative Distribution Function")


###################################################
### code chunk number 9: probhat.Rnw:201-202
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (dFht, main="Quantile Function")


###################################################
### code chunk number 10: probhat.Rnw:211-214
###################################################
dfh (10)
dfh (10, freq=TRUE)
dFht (c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 11: probhat.Rnw:232-235
###################################################
cfh <- pdfuv.cks (height)
cFh <- cdfuv.cks (height)
cFht <- qfuv.cks (height)


###################################################
### code chunk number 12: probhat.Rnw:239-240
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh, TRUE, main="Probability Density Function")


###################################################
### code chunk number 13: probhat.Rnw:242-243
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFh, TRUE, main="Cumulative Distribution Function")


###################################################
### code chunk number 14: probhat.Rnw:245-246
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFht, main="Quantile Function")


###################################################
### code chunk number 15: probhat.Rnw:253-255
###################################################
cfh (22)
cFht (c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 16: probhat.Rnw:266-269
###################################################
#note, depth is third variable
#(will fit a semi-bounded model)
head (quakes2, 2)


###################################################
### code chunk number 17: probhat.Rnw:272-275
###################################################
cfh2 <- pdfmv.cks (quakes2 [,1:2], smoothness = c (0.35, 1) )
cfh3 <- pdfmv.cks (quakes2 [,1:3], smoothness = c (0.35, 1, 1),
    a = c (-Inf, -Inf, 0) )


###################################################
### code chunk number 18: probhat.Rnw:279-281
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2,, TRUE,
    main="Bivariate PDF, 2D")


###################################################
### code chunk number 19: probhat.Rnw:283-284
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2, TRUE, main="Bivariate PDF, 3D")


###################################################
### code chunk number 20: probhat.Rnw:286-288
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh3, main="Trivariate PDF",
    nslides=8, zlim = c (680, 40) )


###################################################
### code chunk number 21: probhat.Rnw:296-298
###################################################
cfh2 (c (180, -20) )
cfh3 (c (180, -20, 300) )


###################################################
### code chunk number 22: probhat.Rnw:314-315
###################################################
conditions <- c (long=180, lat=-20)


###################################################
### code chunk number 23: probhat.Rnw:318-323
###################################################
depth.fhc <- pdfc.cks (quakes2 [,-4], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE,
    a = c (-Inf, -Inf, 0) )
mag.fhc <- pdfc.cks (quakes2 [,-3], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE)


###################################################
### code chunk number 24: probhat.Rnw:327-328
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.fhc, main="conditional distribution of depth\n(long=180, lat=-20)")


###################################################
### code chunk number 25: probhat.Rnw:330-331
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (mag.fhc, main="conditional distribution of mag\n(long=180, lat=-20)")


###################################################
### code chunk number 26: probhat.Rnw:344-350
###################################################
depth.mag.fhc <- pdfmvc.cks (quakes2, smoothness = c (0.35, 1, 1, 1),
    conditions = c (long=180, lat=-20), preserve.range=TRUE,
    a = c (-Inf, -Inf, 0, -Inf) )
lat.long.fhc <- pdfmvc.cks (quakes2 [,-4], smoothness = c (0.35, 1, 1),
    conditions = c (depth=168), preserve.range=TRUE,
    a = c (-Inf, -Inf, 0) )


###################################################
### code chunk number 27: probhat.Rnw:354-356
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.mag.fhc,
    main="conditional distribution of depth and mag\n(long=180, lat=-20)")


###################################################
### code chunk number 28: probhat.Rnw:358-360
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (lat.long.fhc,
    main="conditional distribution of lat and long\n(depth=168)")


###################################################
### code chunk number 29: probhat.Rnw:374-375
###################################################
gfh <- pmfuv.cat (crime.type, n.arrests)


###################################################
### code chunk number 30: probhat.Rnw:379-380
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, main="Probability Mass Function")


###################################################
### code chunk number 31: probhat.Rnw:382-383
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, freq=TRUE, main="same as above\nbut with frequencies")


###################################################
### code chunk number 32: probhat.Rnw:390-391
###################################################
levels (crime.type$crime.type)


###################################################
### code chunk number 33: probhat.Rnw:394-396
###################################################
gfh (1)
gfh ("Assault")


###################################################
### code chunk number 34: probhat.Rnw:399-400
###################################################
gfh ("Assault", freq=TRUE)


###################################################
### code chunk number 35: probhat.Rnw:416-418
###################################################
eFh <- cdf.el (height)
eFht <- qf.el (height)


###################################################
### code chunk number 36: probhat.Rnw:422-423
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFh, main="Cumulative Distribution Function")


###################################################
### code chunk number 37: probhat.Rnw:425-426
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFht, main="Quantile Function")


###################################################
### code chunk number 38: probhat.Rnw:453-457
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1a.gmix <- pmfc.gmix (species, sepal.length,
    conditions = c (sepal.length=5.5) )
plot (fh1a.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=5.5)")


###################################################
### code chunk number 39: probhat.Rnw:459-463
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1b.gmix <- pmfc.gmix (species, sepal.length,
    conditions = c (sepal.length=6.5) )
plot (fh1b.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=6.5)")


###################################################
### code chunk number 40: probhat.Rnw:471-479
###################################################
getOption("SweaveHooks")[["fig"]]()
fh2.gmix <- pmfc.gmix (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )
plot (fh2.gmix,
    main = paste (
        "Conditional Distribution of Iris Species",
        "(sepal.length=6, sepal.width=3)",
        sep="\n")
    )


###################################################
### code chunk number 41: probhat.Rnw:488-490
###################################################
ph.mode (fh2.gmix)
ph.mode (fh2.gmix, level.names=TRUE)


###################################################
### code chunk number 42: probhat.Rnw:508-509
###################################################
fh.gset <- pdfuv.gset.cks (species, sepal.length)


###################################################
### code chunk number 43: probhat.Rnw:513-514
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh.gset, main="Density Estimates of Sepal Length\n(grouped by species)")


###################################################
### code chunk number 44: probhat.Rnw:520-521
###################################################
Fht.mset = qfuv.mset.el (trees2)


###################################################
### code chunk number 45: probhat.Rnw:525-526
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fht.mset, nr=2, nc=2)


###################################################
### code chunk number 46: probhat.Rnw:531-532
###################################################
pdist (fh.gset)


###################################################
### code chunk number 47: probhat.Rnw:548-551
###################################################
#multivariate cdf
#note pwith function can use pdf or cdf, but only nonconditional
cFh3 <- cdfmv.cks (trees2)


###################################################
### code chunk number 48: probhat.Rnw:554-560
###################################################
xlim <- matrix (c (
    22, 24,    #height in 22 to 24
    28, 38,    #girth  in 28 to 38
    0.55, 1.05 #volume in 0.55 to 1.05
    ),, 2, byrow=TRUE, dimnames = list (colnames (trees2), c ("a", "b") ) )
xlim


###################################################
### code chunk number 49: probhat.Rnw:563-566
###################################################
#multivariate probability
probmv (cFh3, xlim [,1], xlim [,2])
pwith (cFh3, xlim=xlim)


###################################################
### code chunk number 50: probhat.Rnw:585-587
###################################################
chFht <- chqf.cks (trees2)
synthetic.data <- rng (chFht, 31)


###################################################
### code chunk number 51: probhat.Rnw:590-598
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
### code chunk number 52: probhat.Rnw:602-604
###################################################
getOption("SweaveHooks")[["fig"]]()
#original data
plot.trees.data (trees2, "original data")


###################################################
### code chunk number 53: probhat.Rnw:607-609
###################################################
getOption("SweaveHooks")[["fig"]]()
#synthetic data
plot.trees.data (synthetic.data, "synthetic data")


###################################################
### code chunk number 54: probhat.Rnw:646-650
###################################################
ph.mean (cFh)
ph.var (cFh)
ph.skewness (cFh)
ph.kurtosis (cFh)


###################################################
### code chunk number 55: probhat.Rnw:662-664
###################################################
quartiles (cFht)
deciles (cFht)


###################################################
### code chunk number 56: probhat.Rnw:667-669
###################################################
ntiles (8, cFht)
ntiles (8, cFht, prob=TRUE)


###################################################
### code chunk number 57: probhat.Rnw:681-683
###################################################
ph.median (cFht)
ph.quantile (cFht, c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 58: probhat.Rnw:686-688
###################################################
#inter-quartile range
iqr (cFht)


###################################################
### code chunk number 59: probhat.Rnw:691-693
###################################################
#inter-quantile ranges
iqr (cFht, 2/3)


###################################################
### code chunk number 60: probhat.Rnw:701-703
###################################################
cFht (0.5)
ph.median (cFht)


###################################################
### code chunk number 61: probhat.Rnw:707-708
###################################################
ph.median (height)


###################################################
### code chunk number 62: probhat.Rnw:720-722
###################################################
ph.mode(cfh)
ph.mode(cfh, TRUE)


###################################################
### code chunk number 63: probhat.Rnw:729-738
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
### code chunk number 64: probhat.Rnw:741-746
###################################################
strs <- c (c ("mean", "median", "mode") )
x <- height.summary [strs]
y <- c (0.06, 0.1, 0.14)
colors <- c ("black", "blue", "darkgreen")
adjv <- c (1.25, 0.5, -0.25)


###################################################
### code chunk number 65: probhat.Rnw:750-754
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh)
abline (v=x, col=colors)
for (i in 1:3)
    text (x [i], y [i], strs [i], adj = adjv [i], col = colors [i])


###################################################
### code chunk number 66: probhat.Rnw:758-759
###################################################
height.summary


###################################################
### code chunk number 67: probhat.Rnw:1079-1080
###################################################
membership <- FKM.gk (unemployment, k=3, seed=2)$U


###################################################
### code chunk number 68: probhat.Rnw:1085-1087
###################################################
w <- membership [,1]
w <- w / sum (w)


###################################################
### code chunk number 69: probhat.Rnw:1092-1093
###################################################
wfh.1 <- pdfmv.cks (unemployment, w=w)


###################################################
### code chunk number 70: probhat.Rnw:1097-1098
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (wfh.1,, TRUE)


###################################################
### code chunk number 71: probhat.Rnw:1101-1103
###################################################
getOption("SweaveHooks")[["fig"]]()
k = 1 - w / max (w)
plot (unemployment, pch=16, col=rgb (k, k, k) )


###################################################
### code chunk number 72: probhat.Rnw:1108-1112
###################################################
w <- membership [,2]
wfh.2 = pdfmv.cks (unemployment, w = w / sum (w) )
w <- membership [,3]
wfh.3 = pdfmv.cks (unemployment, w = w / sum (w) )


###################################################
### code chunk number 73: probhat.Rnw:1116-1119 (eval = FALSE)
###################################################
## plot (wfh.1, main="cluster 1")
## plot (wfh.2, main="cluster 2")
## plot (wfh.3, main="cluster 3")


###################################################
### code chunk number 74: probhat.Rnw:1123-1129
###################################################
getOption("SweaveHooks")[["fig"]]()
p0 <- par (mfrow = c (2, 2) )
plot (wfh.1, main="cluster 1")
plot.new ()
plot (wfh.2, main="cluster 2")
plot (wfh.3, main="cluster 3")
par (p0)


###################################################
### code chunk number 75: probhat.Rnw:1137-1138
###################################################
prep.ph.data (eval=FALSE, echo=TRUE)


###################################################
### code chunk number 76: probhat.Rnw:1146-1150
###################################################
ht <- function (object)
{   print (head (object, 2) )
    print (tail (object, 2) )
}


###################################################
### code chunk number 77: probhat.Rnw:1153-1154
###################################################
ht (cbind (traffic.bins, freq=traffic.freq) )


###################################################
### code chunk number 78: probhat.Rnw:1157-1163
###################################################
ht (trees2)
ht (quakes2)
ht (crimes)
ht (data.frame (crime.type=crime.type$crime.type, n.arrests) )
ht (data.frame (species, sepal.length, sepal.width) )
ht (unemployment)


