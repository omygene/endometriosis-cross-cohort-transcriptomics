## ----style, eval=TRUE, echo=FALSE, results="asis"-----------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE, error=FALSE, tidy=FALSE) # turn off verbosity

## -----------------------------------------------------------------------------
knitr::include_graphics("pcaplot.png")
knitr::include_graphics("arraymetadatatable.png")

## ----options, include=FALSE---------------------------------------------------
##options(error = recover, warn = 2)
options(bitmapType = "cairo")

## ----DataLoading--------------------------------------------------------------
library("ALLMLL")
data("MLL.A")

## ----AffyBatchQM, results="hide"----------------------------------------------
library("arrayQualityMetrics")
arrayQualityMetrics(expressionset = MLL.A[, 1:5],
                    outdir = "Report_for_MLL_A",
                    force = TRUE,
                    do.logtransform = TRUE)

## ----Normalisation, results="hide"--------------------------------------------
nMLL = rma(MLL.A)

## ----ExpressionSet------------------------------------------------------------
arrayQualityMetrics(expressionset = nMLL,
                    outdir = "Report_for_nMLL",
                    force = TRUE)

## ----NChannelSet1, results="hide"---------------------------------------------
library("vsn")
library("CCl4")
data("CCl4")
nCCl4 = justvsn(CCl4, subsample = 15000)
arrayQualityMetrics(expressionset = nCCl4,
                    outdir = "Report_for_nCCl4",
                    force = TRUE)

## ----intgroup1----------------------------------------------------------------
pData(nMLL)$condition = rep(letters[1:4], times = 5)
pData(nMLL)$batch = rep(paste(1:4), each = 5)

## ----intgroup2----------------------------------------------------------------
arrayQualityMetrics(expressionset = nMLL,
                    outdir = "Report_for_nMLL_with_factors",
                    force = TRUE,
                    intgroup = c("condition", "batch"))

## ----XYcoordinates------------------------------------------------------------
featureData(nCCl4)$X = featureData(nCCl4)$Row
featureData(nCCl4)$Y = featureData(nCCl4)$Column

## ----hasTarget----------------------------------------------------------------
featureData(nCCl4)$hasTarget = (regexpr("^NM", featureData(nCCl4)$Name) > 0)
table(featureData(nCCl4)$hasTarget)

## ----pData--------------------------------------------------------------------
pd = pData(CCl4)
rownames(pd) = NULL
pd

## ----RIN----------------------------------------------------------------------
RIN = with(pd, ifelse( Cy3=="CCl4", RIN.Cy3, RIN.Cy5))
fRIN = factor(RIN)
levels(fRIN) = c("poor", "medium", "good")
pData(nCCl4)$"RNA-integrity" = fRIN

## ----NChannelSet2-------------------------------------------------------------
arrayQualityMetrics(expressionset = nCCl4,
                    outdir = "Report_for_nCCl4_with_RIN",
                    force = TRUE,
                    intgroup = "RNA-integrity")

## ----pkgs, echo=FALSE, results="asis"-----------------------------------------
sessionInfo()

