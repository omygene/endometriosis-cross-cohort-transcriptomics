## ----style, eval=TRUE, echo=FALSE, results="asis"-----------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE, error=FALSE, tidy=FALSE) # turn off verbosity

## ----dummy1, include=FALSE----------------------------------------------------
options(bitmapType = "cairo")

## ----loading, results="hide"--------------------------------------------------
library("arrayQualityMetrics")
library("ALLMLL")
data("MLL.A")

## ----DataPreparation----------------------------------------------------------
preparedData = prepdata(expressionset = MLL.A, 
                        intgroup = c(), 
                        do.logtransform = TRUE)

## ----boxplot, eval=FALSE------------------------------------------------------
# ?aqm.boxplot

## ----metrics------------------------------------------------------------------
bo = aqm.boxplot(preparedData)
de = aqm.density(preparedData)
qm = list("Boxplot" = bo, "Density" = de)

## ----booutliers---------------------------------------------------------------
bo@outliers

## ----shortReport--------------------------------------------------------------
outdir = tempdir()
aqm.writereport(modules = qm, reporttitle = "My example", outdir = outdir, 
                arrayTable = pData(MLL.A))
outdir

## ----pkgs, echo=FALSE---------------------------------------------------------
sessionInfo()

