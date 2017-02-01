#' Evaluation of RIVER
#'
#' \code{evaRIVER} trains RIVER by holding out a list of individual and gene
#'         pairs having same rare variants for evaluation, computes test
#'         posterior probabilities of FR for 1st individual, and compares
#'         them with outlier status of 2nd individual from the list.
#'
#' @param dataInput An object of ExpressionSet class which contains input data
#'         required for all functions in RIVERpkg including genomic features,
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
#'         package = "RIVERpkg"), ZscoreThrd=1.5)
#' evaROC <- evaRIVER(dataInput, verbose=TRUE)
#'
#' @export

evaRIVER <- function(dataInput, pseudoc=50,
                     theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Extract required data for evaRIVER
  GAll = t(exprs(dataInput)) # all genomic features (G)
  EAll = as.numeric(unlist(dataInput$Outlier))-1 # all outlier status (E)
  Gtrng = t(exprs(dataInput[,is.na(dataInput$N2pair)])) # G for training models
  Etrng = as.numeric(unlist(dataInput$Outlier
                            [is.na(dataInput$N2pair)]))-1 # E for training models
  # G for test
  Gtest = t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                  exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))
  # E for test (1st and then 2nd individuals from N2 pairs)
  Etest1 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  # E for test (2nd and then 1st individuals from N2 pairs)
  Etest2 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1

  ## Standardization
  mean_g = apply(GAll, 2, mean)
  sd_g = apply(GAll,2,sd)
  GAll <- scale(GAll, center=mean_g, scale=sd_g)
  Gtrng <- scale(Gtrng, center=mean_g, scale=sd_g)

  ## Search a best lambda from a multivariate logistic regression
  ##         with outlier status with 10 cross-validation
  ## GAM (genomeic annotation model)
  cv.ll = cv.glmnet(Gtrng, as.vector(Etrng), lambda=costs,
                    family="binomial", alpha=0, nfolds=10)
  if (verbose) cat(' *** best lambda = ',cv.ll$lambda.min,' *** \n\n', sep='')

  ## Compute a P(FR | G) for all data
  pp_test = predict(cv.ll, Gtest, s="lambda.min", type="response") # >>

  ## Train RIVER on training data
  em.res <- integratedEM(Gtrng, Etrng, cv.ll$lambda.min,
                         cv.ll$glmnet.fit, pseudoc, theta_init, costs, verbose)

  # ## Generate G data for test data (Revised)
  Gtest <- scale(Gtest, center=mean_g, scale=sd_g)

  ## Compute P(FR | G, E)
  dup.post = testPosteriors(Gtest, Etest1, em.res)

  ## Check performance of models with N2 pairs
  RIVER.roc = roc(Etest2, dup.post$posterior[,2]) # RIVER
  GAM.roc = roc(Etest2, as.numeric(pp_test)) # GAM

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
