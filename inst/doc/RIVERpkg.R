## ----'ultraQuick'--------------------------------------------------------
library("RIVERpkg")

## ------------------------------------------------------------------------
data(simulated_features) # G: features from genomic annotations
data(simulated_outliers) # E: outlier status from expression data

## ------------------------------------------------------------------------
head(simulated_features)

## ------------------------------------------------------------------------
head(simulated_outliers)

## ------------------------------------------------------------------------
rocSTAT <- evaRIVER(simulated_features, simulated_outliers)

## ------------------------------------------------------------------------
summary(rocSTAT)

## ------------------------------------------------------------------------
plotAUC(rocSTAT)

## ------------------------------------------------------------------------
outRIVER <- appRIVER(simulated_features, simulated_outliers)

## ------------------------------------------------------------------------
summary(outRIVER)

## ------------------------------------------------------------------------
plotPosteriors(outRIVER, simulated_outliers)

## ------------------------------------------------------------------------
rocSTAT <- evaRIVER(simulated_features, simulated_outliers, pseudoc=50, 
                    theta_init=matrix(c(.99, .01, .3, .7), nrow=2), 
                    costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=TRUE)

## ------------------------------------------------------------------------
outRIVER <- appRIVER(simulated_features, simulated_outliers, pseudoc=50, 
                     theta_init=matrix(c(.99, .01, .3, .7), nrow=2), 
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=TRUE)

## ------------------------------------------------------------------------
print(outRIVER$fitRIVER$beta)

## ------------------------------------------------------------------------
print(outRIVER$fitRIVER$theta)

## ---- eval = FALSE-------------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("RIVERpkg")

## ---------------------------------------------------------------------------------------------------------------------
## Session info
library('devtools')
options(width = 120)
session_info()

