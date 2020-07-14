### R code from vignette source 'probhat.Rnw'

###################################################
### code chunk number 1: probhat.Rnw:35-39
###################################################
options(continue=" ", width=80)
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 2.6, 1.6), cex=0.7, cex.main=1)))
set.seed (1)


###################################################
### code chunk number 2: probhat.Rnw:153-157
###################################################
library (probhat)       #required
library (vectools)      #optional, used in appendix
library (fclust)        #optional, used in appendix
library (scatterplot3d) #optional


###################################################
### code chunk number 3: probhat.Rnw:164-165
###################################################
set.ph.options (rendering.style="e", theme="green")


###################################################
### code chunk number 4: probhat.Rnw:170-171
###################################################
ph.data.prep ()


###################################################
### code chunk number 5: probhat.Rnw:185-188
###################################################
dfh <- pmfuv.dks (traffic.bins, traffic.freq, lower=0)
dFh <- cdfuv.dks (traffic.bins, traffic.freq, lower=0)
dFht <- qfuv.dks (traffic.bins, traffic.freq, lower=0)


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
### code chunk number 11: probhat.Rnw:231-234
###################################################
cfh <- pdfuv.cks (height)
cFh <- cdfuv.cks (height)
cFht <- qfuv.cks (height)


###################################################
### code chunk number 12: probhat.Rnw:238-239
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh, TRUE, main="Probability Density Function")


###################################################
### code chunk number 13: probhat.Rnw:241-242
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFh, TRUE, main="Cumulative Distribution Function")


###################################################
### code chunk number 14: probhat.Rnw:244-245
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cFht, main="Quantile Function")


###################################################
### code chunk number 15: probhat.Rnw:252-254
###################################################
cfh (22)
cFht (c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 16: probhat.Rnw:265-267
###################################################
cfh2 <- pdfmv.cks (quakes [,1:2], smoothness = c (0.35, 1) )
cfh3 <- pdfmv.cks (quakes [,1:3], smoothness = c (0.35, 1, 1) )


###################################################
### code chunk number 17: probhat.Rnw:271-273
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2,, TRUE,
    main="Bivariate PDF, 2D")


###################################################
### code chunk number 18: probhat.Rnw:275-276
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh2, TRUE, main="Bivariate PDF, 3D")


###################################################
### code chunk number 19: probhat.Rnw:278-280
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh3, main="Trivariate PDF",
    nslides=8, zlim = c (680, 40) )


###################################################
### code chunk number 20: probhat.Rnw:288-290
###################################################
cfh2 (c (180, -20) )
cfh3 (c (180, -20, 300) )


###################################################
### code chunk number 21: probhat.Rnw:306-307
###################################################
conditions <- c (long=180, lat=-20)


###################################################
### code chunk number 22: probhat.Rnw:310-314
###################################################
depth.fhc <- pdfc.cks (quakes [,-4], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE)
mag.fhc <- pdfc.cks (quakes [,-3], smoothness = c (0.35, 1, 1),
    conditions=conditions, preserve.range=TRUE)


###################################################
### code chunk number 23: probhat.Rnw:318-319
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.fhc, main="conditional distribution of depth\n(long=180, lat=-20)")


###################################################
### code chunk number 24: probhat.Rnw:321-322
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (mag.fhc, main="conditional distribution of mag\n(long=180, lat=-20)")


###################################################
### code chunk number 25: probhat.Rnw:335-339
###################################################
depth.mag.fhc <- pdfmvc.cks (quakes, smoothness = c (0.35, 1, 1, 1),
    conditions = c (long=180, lat=-20), preserve.range=TRUE)
lat.long.fhc <- pdfmvc.cks (quakes [,-4], smoothness = c (0.35, 1, 1),
    conditions = c (depth=168), preserve.range=TRUE)


###################################################
### code chunk number 26: probhat.Rnw:343-345
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (depth.mag.fhc,
    main="conditional distribution of depth and mag\n(long=180, lat=-20)")


###################################################
### code chunk number 27: probhat.Rnw:347-349
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (lat.long.fhc,
    main="conditional distribution of lat and long\n(depth=168)")


###################################################
### code chunk number 28: probhat.Rnw:363-364
###################################################
gfh <- pmfuv.cat (crime.type, n.arrests)


###################################################
### code chunk number 29: probhat.Rnw:368-369
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, main="Probability Mass Function")


###################################################
### code chunk number 30: probhat.Rnw:371-372
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (gfh, freq=TRUE, main="same as above\nbut with frequencies")


###################################################
### code chunk number 31: probhat.Rnw:379-380
###################################################
levels (crime.type$crime.type)


###################################################
### code chunk number 32: probhat.Rnw:383-385
###################################################
gfh (1)
gfh ("Assault")


###################################################
### code chunk number 33: probhat.Rnw:388-389
###################################################
gfh ("Assault", freq=TRUE)


###################################################
### code chunk number 34: probhat.Rnw:405-407
###################################################
eFh <- cdf.el (height)
eFht <- qf.el (height)


###################################################
### code chunk number 35: probhat.Rnw:411-412
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFh, main="Cumulative Distribution Function")


###################################################
### code chunk number 36: probhat.Rnw:414-415
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (eFht, main="Quantile Function")


###################################################
### code chunk number 37: probhat.Rnw:442-446
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1a.gmix <- pmfc.gmixp (species, sepal.length,
    conditions = c (sepal.length=5.5) )
plot (fh1a.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=5.5)")


###################################################
### code chunk number 38: probhat.Rnw:448-452
###################################################
getOption("SweaveHooks")[["fig"]]()
fh1b.gmix <- pmfc.gmixp (species, sepal.length,
    conditions = c (sepal.length=6.5) )
plot (fh1b.gmix,
    main="Conditional Distribution of Iris Species\n(sepal.length=6.5)")


###################################################
### code chunk number 39: probhat.Rnw:460-468
###################################################
getOption("SweaveHooks")[["fig"]]()
fh2.gmix <- pmfc.gmixp (species, cbind (sepal.length, sepal.width),
    conditions = c (sepal.length=6, sepal.width=3) )
plot (fh2.gmix,
    main = paste (
        "Conditional Distribution of Iris Species",
        "(sepal.length=6, sepal.width=3)",
        sep="\n")
    )


###################################################
### code chunk number 40: probhat.Rnw:477-479
###################################################
ph.mode (fh2.gmix)
gmode (fh2.gmix)


###################################################
### code chunk number 41: probhat.Rnw:495-496
###################################################
fh.gset <- pdfuv.gset.cks (sepal.length, group.by=species)


###################################################
### code chunk number 42: probhat.Rnw:500-501
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh.gset, main="Density Estimates of Sepal Length\n(grouped by species)")


###################################################
### code chunk number 43: probhat.Rnw:507-508
###################################################
Fht.mset = qfuv.mset.el (trees)


###################################################
### code chunk number 44: probhat.Rnw:512-513
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fht.mset, nr=2, nc=2)


###################################################
### code chunk number 45: probhat.Rnw:525-527
###################################################
#multivariate cdf
cFh3 <- cdfmv.cks (trees)


###################################################
### code chunk number 46: probhat.Rnw:530-536
###################################################
q <- matrix (c (
    22, 24,    #height in 22 to 24
    28, 38,    #girth  in 28 to 38
    0.55, 1.05 #volume in 0.55 to 1.05
    ),, 2, byrow=TRUE, dimnames = list (colnames (trees), c ("a", "b") ) )
q


###################################################
### code chunk number 47: probhat.Rnw:539-541
###################################################
#multivariate probability
probmv (cFh3, q [,1], q [,2])


###################################################
### code chunk number 48: probhat.Rnw:560-562
###################################################
chFht <- chqf.cks (trees)
synthetic.data <- rng (chFht, 31)


###################################################
### code chunk number 49: probhat.Rnw:565-573
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
### code chunk number 50: probhat.Rnw:577-579
###################################################
getOption("SweaveHooks")[["fig"]]()
#original data
plot.trees.data (trees, "original data")


###################################################
### code chunk number 51: probhat.Rnw:582-584
###################################################
getOption("SweaveHooks")[["fig"]]()
#synthetic data
plot.trees.data (synthetic.data, "synthetic data")


###################################################
### code chunk number 52: probhat.Rnw:621-625
###################################################
ph.mean (cFh)
ph.var (cFh)
ph.skewness (cFh)
ph.kurtosis (cFh)


###################################################
### code chunk number 53: probhat.Rnw:637-639
###################################################
quartiles (cFht)
deciles (cFht)


###################################################
### code chunk number 54: probhat.Rnw:642-644
###################################################
ntiles (cFht, 8)
ntiles (cFht, 8, rank=FALSE)


###################################################
### code chunk number 55: probhat.Rnw:647-650
###################################################
quartiles (cFht)
quartiles (eFht)
quartiles (height)


###################################################
### code chunk number 56: probhat.Rnw:662-664
###################################################
ph.median (cFht)
ph.quantile (cFht, c (0.25, 0.5, 0.75) )


###################################################
### code chunk number 57: probhat.Rnw:667-669
###################################################
#inter-quartile range
iqr (cFht)


###################################################
### code chunk number 58: probhat.Rnw:672-674
###################################################
#inter-quantile ranges
iqr (cFht, 2/3)


###################################################
### code chunk number 59: probhat.Rnw:682-684
###################################################
cFht (0.5)
ph.median (cFht)


###################################################
### code chunk number 60: probhat.Rnw:688-689
###################################################
ph.median (height)


###################################################
### code chunk number 61: probhat.Rnw:701-703
###################################################
ph.mode(cfh)
ph.mode(cfh, TRUE)


###################################################
### code chunk number 62: probhat.Rnw:710-719
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
### code chunk number 63: probhat.Rnw:722-727
###################################################
strs <- c (c ("mean", "median", "mode") )
x <- height.summary [strs]
y <- c (0.06, 0.1, 0.14)
colors <- c ("black", "blue", "darkgreen")
adjv <- c (1.25, 0.5, -0.25)


###################################################
### code chunk number 64: probhat.Rnw:731-735
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh)
abline (v=x, col=colors)
for (i in 1:3)
    text (x [i], y [i], strs [i], adj = adjv [i], col = colors [i])


###################################################
### code chunk number 65: probhat.Rnw:739-740
###################################################
height.summary


###################################################
### code chunk number 66: probhat.Rnw:1146-1147
###################################################
membership <- FKM.gk (unemployment, k=3, seed=2)$U


###################################################
### code chunk number 67: probhat.Rnw:1152-1154
###################################################
w <- membership [,1]
w <- w / sum (w)


###################################################
### code chunk number 68: probhat.Rnw:1159-1160
###################################################
wfh.1 <- pdfmv.cks (unemployment, w=w)


###################################################
### code chunk number 69: probhat.Rnw:1164-1165
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (wfh.1,, TRUE)


###################################################
### code chunk number 70: probhat.Rnw:1168-1170
###################################################
getOption("SweaveHooks")[["fig"]]()
k = 1 - w / max (w)
plot (unemployment, pch=16, col=rgb (k, k, k) )


###################################################
### code chunk number 71: probhat.Rnw:1175-1179
###################################################
w <- membership [,2]
wfh.2 = pdfmv.cks (unemployment, w = w / sum (w) )
w <- membership [,3]
wfh.3 = pdfmv.cks (unemployment, w = w / sum (w) )


###################################################
### code chunk number 72: probhat.Rnw:1183-1186 (eval = FALSE)
###################################################
## plot (wfh.1, main="cluster 1")
## plot (wfh.2, main="cluster 2")
## plot (wfh.3, main="cluster 3")


###################################################
### code chunk number 73: probhat.Rnw:1190-1196
###################################################
getOption("SweaveHooks")[["fig"]]()
p0 <- par (mfrow = c (2, 2) )
plot (wfh.1, main="cluster 1")
plot.new ()
plot (wfh.2, main="cluster 2")
plot (wfh.3, main="cluster 3")
par (p0)


###################################################
### code chunk number 74: probhat.Rnw:1204-1205
###################################################
ph.data.prep (eval=FALSE, echo=TRUE)


###################################################
### code chunk number 75: probhat.Rnw:1213-1214
###################################################
headt (cbind (traffic.bins, freq=traffic.freq) )


###################################################
### code chunk number 76: probhat.Rnw:1217-1223
###################################################
headt (trees)
headt (quakes)
headt (crimes)
headt (data.frame (crime.type=crime.type$crime.type, n.arrests) )
headt (data.frame (species, sepal.length, sepal.width) )
headt (unemployment)


