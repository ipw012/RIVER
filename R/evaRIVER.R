#' Evaluation of RIVER
#'
#' \code{evaRIVER} trains RIVER by holding out a list of individual and gene
#'         pairs having same rare variants for evaluation, computes test
#'         posterior probabilities of FR for 1st individual, and compares
#'         them with outlier status of 2nd individual from the list.
#'
#' @param dataInput An object of ExpressionSet class which contains input data
#'         required for all functions in RIVER including genomic features,
#'         outlier status, and N2 pairs.
#' @param pseudoc Pseudo count.
#' @param theta_init Initial values of theta.
#' @param costs Candidate penalty parameter values for L2-regularized logistic
#'         regression.
#' @param verbose Logical option for showing extra information on progress.
#'
#' @return A list which contains two AUC values from RIVER and GAM, computed
#'         specificities and sensitivities from two models, and P-value of
#'         comparing the two AUC values.
#'
#' @section Warning: A vector of candidate penalty values make \code{glmnet}
#'         faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{predict}},
#'         \code{\link{integratedEM}}, \code{\link{testPosteriors}},
#'         \code{\link{getData}}, \code{\link[Biobase]{exprs}}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' evaROC <- evaRIVER(dataInput, verbose=TRUE)
#'
#' @export

evaRIVER <- function(dataInput, pseudoc=50,
                     theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Extract required data for evaRIVER
  FeatAll = t(exprs(dataInput)) # all genomic features (G)
  OutAll = as.numeric(unlist(dataInput$Outlier))-1 # all outlier status (E)
  FeatTrng = t(exprs(dataInput[,is.na(dataInput$N2pair)])) # G for training models
  OutTrng = as.numeric(unlist(dataInput$Outlier
                            [is.na(dataInput$N2pair)]))-1 # E for training models
  # G for test
  FeatTest = t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                  exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))
  # E for test (1st and then 2nd individuals from N2 pairs)
  OutTest1 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  # E for test (2nd and then 1st individuals from N2 pairs)
  OutTest2 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1

  ## Standardization
  meanFeat = apply(FeatAll, 2, mean)
  sdFeat = apply(FeatAll,2,sd)
  FeatAll <- scale(FeatAll, center=meanFeat, scale=sdFeat)
  FeatTrng <- scale(FeatTrng, center=meanFeat, scale=sdFeat)

  ## Search a best lambda from a multivariate logistic regression
  ##         with outlier status with 10 cross-validation
  ## GAM (genomeic annotation model)
  logisticCV = cv.glmnet(FeatTrng, as.vector(OutTrng), lambda=costs,
                    family="binomial", alpha=0, nfolds=10)
  if (verbose) cat(' *** best lambda = ',logisticCV$lambda.min,' *** \n\n', sep='')

  ## Compute a P(FR | G) for all data
  postprobTest = predict(logisticCV, FeatTest, s="lambda.min", type="response") # >>

  ## Train RIVER on training data
  emModel <- integratedEM(FeatTrng, OutTrng, logisticCV$lambda.min,
                         logisticCV$glmnet.fit, pseudoc, theta_init, costs, verbose)

  # ## Generate G data for test data (Revised)
  FeatTest <- scale(FeatTest, center=meanFeat, scale=sdFeat)

  ## Compute P(FR | G, E)
  dup.post = testPosteriors(FeatTest, OutTest1, emModel)

  ## Check performance of models with N2 pairs
  RIVER.roc = roc(OutTest2, dup.post$posterior[,2]) # RIVER
  GAM.roc = roc(OutTest2, as.numeric(postprobTest)) # GAM

  if (verbose) {
    cat('*** AUC (GAM - genomic annotation model): ',round(GAM.roc$auc,3),'\n')
    cat('    AUC (RIVER): ',round(RIVER.roc$auc,3),'\n')
    cat('    P-value: ',format.pval(roc.test(RIVER.roc, GAM.roc)$p.value,digits=2,
                                    eps=0.001),'***\n')
    cat('\n')
  }

  evaROC=list(RIVER_sens=RIVER.roc$sensitivities, RIVER_spec=RIVER.roc$specificities,
       RIVER_auc=RIVER.roc$auc[1], GAM_sens=GAM.roc$sensitivities,
       GAM_spec=GAM.roc$specificities, GAM_auc=GAM.roc$auc[1],
       pvalue=roc.test(RIVER.roc, GAM.roc)$p.value)
  class(evaROC) = "eval"
  return(evaROC)
}
