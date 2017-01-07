#' Evaluation of RIVER
#'
#' \code{evaRIVER} trains RIVER by holding out a list of individual and gene
#'         pairs having same rare variants for evaluation, computes test
#'         posterior probabilities of FR for 1st individual, and compares
#'         them with outlier status of 2nd individual from the list.
#'
#' @param G Genomic features in the data frame which contains subject IDs,
#'         gene symbols, and P number of genomic fetures of each instance
#'         per row.
#' @param E Expression outliers in the data frame which contains subject IDs,
#'         gene symbols, and outlier status of each instance per row.
#' @param pseudoc Pseudo count.
#' @param theta_init Initial values of theta.
#' @param costs Candidate penalty parameter values for L2-regularized logistic
#'         regression.
#' @param verbose Logical option for showing extra information on progress.
#'
#' @return An object of class RIVER which include two AUC values from
#'         RIVER and GAM and corresponding P-value of comparing the two AUC
#'         values.
#'
#' @section Warning: A vector of candidate penalty values make \code{glmnet}
#'         faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{predict}},
#'         \code{\link{integratedEM}}, and \code{\link{testPosteriors}}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' rocSTAT <- evaRIVER(G, E, verbose=TRUE)
#'
#' @export

evaRIVER <- function(G, E, pseudoc=50, theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Search N2 pairs
  ## (a list of two individuals having same rare SNVs within 10Kb of TSS)
  grp.cols <- names(G)[2:ncol(G)]
  dots <- lapply(grp.cols, as.symbol)
  K = data.frame(G %>% mutate(rnum = row_number()) %>% group_by_(.dots = dots)
                 %>% arrange() %>% mutate(number=n()) %>% filter(number == 2))
  dups <- matrix(K$rnum, length(K$rnum)/2, 2, byrow = TRUE)

  E_disc = as.numeric(E[,"Outlier"])

  basic_data = cbind(E_disc, G)
  rp = sample.int(dim(basic_data)[1])

  ## Split data into trainig and test dataset
  train_inds = setdiff(rp, union(dups[,1], dups[,2]))
  test_inds = union(dups[,1], dups[,2])

  g_all <- scale(as.matrix(G[,3:ncol(G)]))

  ## Generate training data by holding out N2 pairs
  g_trng <- scale(as.matrix(G[train_inds,3:ncol(G)]),
                  center=colMeans(g_all), scale=apply(g_all,2,sd))

  ## Search a best lambda from a multivariate logistic regression
  ##         with outlier status with 10 cross-validation

  ## GAM (genomeic annotation model)
  cv.ll = cv.glmnet(g_trng, as.vector(E_disc[train_inds]), lambda=costs,
                    family="binomial", alpha=0, nfolds=10)
  if (verbose) cat(' *** best lambda = ',cv.ll$lambda.min,' *** \n\n', sep='')

  ## Compute a P(FR | G) for all data
  pp = predict(cv.ll, g_all, s="lambda.min", type="response")

  ## Train RIVER on training data
  em.res <- integratedEM(g_trng, E_disc[train_inds], cv.ll$lambda.min,
                         cv.ll$glmnet.fit, pseudoc, theta_init, costs, verbose)

  ## Generate G data for test data (individual 1 from N2 pairs)
  g_test <- scale(as.matrix(G[dups[,2], 3:ncol(G)]),
                  center=colMeans(g_all), scale=apply(g_all,2,sd))

  ## Compute P(FR | G, E)
  dup.post = testPosteriors(g_test, E_disc[dups[,2]], em.res)

  ## Check performance of models with N2 pairs
  RIVER.roc = roc(basic_data$E_disc[dups[,1]], dup.post$posterior[,2]) # RIVER
  GAM.roc = roc(basic_data$E_disc[dups[,1]], pp[dups[,1]]) # GAM

  if (verbose) {
    cat('*** AUC (GAM - genomic annotation model): ',round(GAM.roc$auc,3),'\n')
    cat('    AUC (RIVER): ',round(RIVER.roc$auc,3),'\n')
    cat('    P-value: ',format.pval(roc.test(RIVER.roc, GAM.roc)$p.value,digits=2,
                                    eps=0.001),'***\n')
    cat('\n')
  }

  list(rocRIVER=RIVER.roc, rocGAM=GAM.roc, pValue=roc.test(RIVER.roc, GAM.roc)$p.value)
}
