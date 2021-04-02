### R code from vignette source 'probhatx_kernels.Rnw'

###################################################
### code chunk number 1: probhatx_kernels.Rnw:19-21
###################################################
options(SweaveHooks=list(fig=function()
par(mar=c(4.1, 4.1, 2.6, 1.6), cex=0.7, cex.main=1)))


###################################################
### code chunk number 2: probhatx_kernels.Rnw:32-34
###################################################
getOption("SweaveHooks")[["fig"]]()
library (probhat)
plot_kernel_array ()


