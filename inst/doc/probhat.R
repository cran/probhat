### R code from vignette source 'probhat.Rnw'

###################################################
### code chunk number 1: probhat.Rnw:35-39
###################################################
options(continue=" ")
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 2.6, 1.6), cex=0.7, cex.main=1)))
set.seed (1)


###################################################
### code chunk number 2: probhat.Rnw:198-201
###################################################
library (probhat)
library (fclust)
library (scatterplot3d)


###################################################
### code chunk number 3: probhat.Rnw:208-209
###################################################
use.ph.theme ("green")


###################################################
### code chunk number 4: probhat.Rnw:214-215
###################################################
data.prep ()


###################################################
### code chunk number 5: probhat.Rnw:227-228
###################################################
fh = pmfuv.dks (traffic.x, traffic.h, bw=23, lower=0)


###################################################
### code chunk number 6: probhat.Rnw:232-233
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh, TRUE)


###################################################
### code chunk number 7: probhat.Rnw:243-244
###################################################
fh (10)


###################################################
### code chunk number 8: probhat.Rnw:249-251
###################################################
Fh = cdfuv.dks (traffic.x, traffic.h, bw=23, lower=0)
Fh.inv = qfuv.dks (traffic.x, traffic.h, bw=23, lower=0)


###################################################
### code chunk number 9: probhat.Rnw:255-256
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh, TRUE)


###################################################
### code chunk number 10: probhat.Rnw:259-260
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh.inv)


###################################################
### code chunk number 11: probhat.Rnw:279-280
###################################################
fh = pdfuv.cks (Height)


###################################################
### code chunk number 12: probhat.Rnw:284-285
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh, TRUE)


###################################################
### code chunk number 13: probhat.Rnw:293-294
###################################################
fh (22)


###################################################
### code chunk number 14: probhat.Rnw:299-301
###################################################
Fh = cdfuv.cks (Height)
Fh.inv = qfuv.cks (Height)


###################################################
### code chunk number 15: probhat.Rnw:305-306
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh, TRUE)


###################################################
### code chunk number 16: probhat.Rnw:309-310
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh.inv)


###################################################
### code chunk number 17: probhat.Rnw:318-320
###################################################
p1 = 0.5
p2 = Fh (Fh.inv (p1) )


###################################################
### code chunk number 18: probhat.Rnw:323-325
###################################################
p1 == p2
abs (p1 - p2)


###################################################
### code chunk number 19: probhat.Rnw:335-336
###################################################
fh = pdfmv.cks (trees [,2:3])


###################################################
### code chunk number 20: probhat.Rnw:339-340
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh, all=TRUE)


###################################################
### code chunk number 21: probhat.Rnw:346-347
###################################################
fh (c (22, 0.8) )


###################################################
### code chunk number 22: probhat.Rnw:351-352
###################################################
Fh = cdfmv.cks (trees [,2:3])


###################################################
### code chunk number 23: probhat.Rnw:356-357
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh, all=TRUE)


###################################################
### code chunk number 24: probhat.Rnw:368-370
###################################################
conditions = c (Girth=30, Height=22)
cfh = pdfc.cks (trees, conditions=conditions)


###################################################
### code chunk number 25: probhat.Rnw:374-375
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh)


###################################################
### code chunk number 26: probhat.Rnw:380-382
###################################################
#density of volume (volume=0.85), given girth=30 and height=22
cfh (0.85)


###################################################
### code chunk number 27: probhat.Rnw:396-398
###################################################
conditions = c (lat=-20, long=180)
cfh = pdfmvc.cks (quakes, conditions=conditions)


###################################################
### code chunk number 28: probhat.Rnw:402-403
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh, xlim = c (0, 800) )


###################################################
### code chunk number 29: probhat.Rnw:405-406
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh, TRUE, xlim = c (0, 800) )


###################################################
### code chunk number 30: probhat.Rnw:416-420
###################################################
cfh.AA = pdfmvc.cks (quakes, conditions = c (lat=-30, long=170) )
cfh.AB = pdfmvc.cks (quakes, conditions = c (lat=-30, long=180) )
cfh.BA = pdfmvc.cks (quakes, conditions = c (lat=-20, long=170) )
cfh.BB = cfh


###################################################
### code chunk number 31: probhat.Rnw:424-428
###################################################
getOption("SweaveHooks")[["fig"]]()
plot_2x2 (cfh.AA, cfh.AB, cfh.BA, cfh.BB,
    "lat=-30, long=170", "lat=-30, long=180",
    "lat=-20, long=170", "lat=-20, long=180",
    xlim = c (0, 800) )


###################################################
### code chunk number 32: probhat.Rnw:434-436
###################################################
cfh.A = pdfmvc.cks (quakes [,-4], conditions = c (depth=100) )
cfh.B = pdfmvc.cks (quakes [,-4], conditions = c (depth=500) )


###################################################
### code chunk number 33: probhat.Rnw:440-441
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh.A, main="depth=100")


###################################################
### code chunk number 34: probhat.Rnw:443-444
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh.B, main="depth=500")


###################################################
### code chunk number 35: probhat.Rnw:457-458
###################################################
fh = pmfuv.cat (region)


###################################################
### code chunk number 36: probhat.Rnw:462-463
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh)


###################################################
### code chunk number 37: probhat.Rnw:468-470
###################################################
fh (1)
fh ("North Central")


###################################################
### code chunk number 38: probhat.Rnw:484-486
###################################################
mean.Sepal.Length = mean (iris.Sepal.Length)
fh = pmfc.cat.cks (iris.Species, iris.Sepal.Length, at=mean.Sepal.Length)


###################################################
### code chunk number 39: probhat.Rnw:490-491
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (fh)


###################################################
### code chunk number 40: probhat.Rnw:504-505
###################################################
Fh = cdf.el (Height)


###################################################
### code chunk number 41: probhat.Rnw:509-510
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh)


###################################################
### code chunk number 42: probhat.Rnw:516-517
###################################################
Fh.inv = qf.el (Height)


###################################################
### code chunk number 43: probhat.Rnw:521-522
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (Fh.inv)


###################################################
### code chunk number 44: probhat.Rnw:550-551
###################################################
cs = categorical.set (pdfuv.cks, iris.Sepal.Length, group.by=iris.Species)


###################################################
### code chunk number 45: probhat.Rnw:555-556
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cs)


###################################################
### code chunk number 46: probhat.Rnw:561-563
###################################################
conditions = cbind (Height = c (20, 24, 28) )
cond.set = conditional.set (pdfc.cks, trees [,2:3], group.by=conditions)


###################################################
### code chunk number 47: probhat.Rnw:567-568
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cond.set)


###################################################
### code chunk number 48: probhat.Rnw:578-579
###################################################
ms = marginal.set (qf.el, trees)


###################################################
### code chunk number 49: probhat.Rnw:583-584
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (ms)


###################################################
### code chunk number 50: probhat.Rnw:596-598
###################################################
#multivariate cdf
Fh = cdfmv.cks (trees)


###################################################
### code chunk number 51: probhat.Rnw:601-604
###################################################
#approximate first and third quartiles
a = c (28, 22, 0.55)
b = c (38, 24, 1.05)


###################################################
### code chunk number 52: probhat.Rnw:607-608
###################################################
cbind (lower=a, upper=b)


###################################################
### code chunk number 53: probhat.Rnw:611-613
###################################################
#multivariate probability
probmv (Fh, a, b)


###################################################
### code chunk number 54: probhat.Rnw:642-644
###################################################
chF.inv = chqf.cks (trees)
synthetic.data = ph.rng (chF.inv, 31)


###################################################
### code chunk number 55: probhat.Rnw:647-655
###################################################
#convenience function
plot.trees.data = function (x, main)
{   Height = x [,"Height"]
    Girth = x [,"Girth"]
    Volume = x [,"Volume"]
    scatterplot3d (Height, Girth, Volume,
        main=main, type="h", angle=112.5, pch=16)
}


###################################################
### code chunk number 56: probhat.Rnw:659-661
###################################################
getOption("SweaveHooks")[["fig"]]()
#original data
plot.trees.data (trees, "original data")


###################################################
### code chunk number 57: probhat.Rnw:664-666
###################################################
getOption("SweaveHooks")[["fig"]]()
#synthetic data
plot.trees.data (synthetic.data, "synthetic data")


###################################################
### code chunk number 58: probhat.Rnw:681-686
###################################################
#conditional distributions
conditions = c (Height=24)
cfh = pdfc.cks (trees [,2:3], conditions=conditions)
cFh = cdfc.cks (trees [,2:3], conditions=conditions)
cFh.inv = qfc.cks (trees [,2:3], conditions=conditions)


###################################################
### code chunk number 59: probhat.Rnw:689-693
###################################################
#mean, median and mode
mean.Volume = ph.mean (cFh)
median.Volume = cFh.inv (0.5)
mode.Volume = ph.mode (cfh)


###################################################
### code chunk number 60: probhat.Rnw:696-698
###################################################
cbind (statistic = c ("mean", "median", "mode"),
    value = c (mean.Volume, median.Volume, mode.Volume) )


###################################################
### code chunk number 61: probhat.Rnw:702-705
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (cfh)
abline (v = c (mean.Volume, mode.Volume) )
abline (v=median.Volume, lty=2)


###################################################
### code chunk number 62: probhat.Rnw:709-713
###################################################
#and just as an example, the variance, skewness and kurtosis...
ph.var (cFh)
ph.skewness (cFh)
ph.kurtosis (cFh)


###################################################
### code chunk number 63: probhat.Rnw:1028-1029
###################################################
data.prep (eval=FALSE, echo=TRUE)


###################################################
### code chunk number 64: probhat.Rnw:1046-1047
###################################################
membership = FKM.gk (unemployment, k=3, seed=2)$U


###################################################
### code chunk number 65: probhat.Rnw:1052-1054
###################################################
w = membership [,1]
w = w / sum (w)


###################################################
### code chunk number 66: probhat.Rnw:1059-1060
###################################################
wfh.1 = pdfmv.cks (unemployment, w=w)


###################################################
### code chunk number 67: probhat.Rnw:1064-1065
###################################################
getOption("SweaveHooks")[["fig"]]()
plot (wfh.1)


###################################################
### code chunk number 68: probhat.Rnw:1068-1070
###################################################
getOption("SweaveHooks")[["fig"]]()
k = 1 - w / max (w)
plot (unemployment, pch=16, col=rgb (k, k, k) )


###################################################
### code chunk number 69: probhat.Rnw:1075-1079
###################################################
w = membership [,2]
wfh.2 = pdfmv.cks (unemployment, w = w / sum (w) )
w = membership [,3]
wfh.3 = pdfmv.cks (unemployment, w = w / sum (w) )


###################################################
### code chunk number 70: probhat.Rnw:1084-1085
###################################################
getOption("SweaveHooks")[["fig"]]()
plot_2x2 (wfh.1,, wfh.2, wfh.3, "cluster 1",, "cluster 2", "cluster 3")


