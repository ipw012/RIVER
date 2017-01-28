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
#' @param P A list of N2 pairs in the data frame which contains pairs of
#'         examples with same rare variants. First two columns contain
#'         SubjectID1 and GeneName1 for one and next two columns contain
#'         SubjectID2 and GeneName2 for another.
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
#' P <- simulated_N2pairs
#' rocSTAT <- evaRIVER(G, E, P, verbose=TRUE)
#'
#' @export

evaRIVER <- function(G, E, P, pseudoc=50, theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Extract indices of N2 pairs (Revised)
  ## (a list of two individuals having same rare SNVs within 10Kb of TSS)
  G2 = data.frame(G %>% mutate(rnum=row_number()))
  temp_N2pairs = P
  colnames(temp_N2pairs) = rep(c("SubjectID","GeneName"),times=2)
  dups = cbind(left_join(temp_N2pairs[,1:2], G2, by=c("SubjectID","GeneName"))[,"rnum"],
                left_join(temp_N2pairs[,3:4], G2, by=c("SubjectID","GeneName"))[,"rnum"])

  E_disc = as.numeric(E[,"Outlier"])

  ## Split data into trainig and test ones
  train_inds = setdiff(sample.int(nrow(G)), union(dups[,1],dups[,2]))
  test_inds = as.integer(rbind(dups[,1],dups[,2]))
  test_inds_rev = as.integer(rbind(dups[,2],dups[,1]))

  ## Standardization
  g_all = as.matrix(G[,3:ncol(G)])
  mean_col_g = apply(g_all,2,mean)
  sd_col_g = apply(g_all,2,sd)
  g_all <- scale(g_all, center=mean_col_g, scale=sd_col_g)

  ## Generate training data by holding out N2 pairs
  g_trng <- scale(as.matrix(G[train_inds,3:ncol(G)]),
                  center=mean_col_g, scale=sd_col_g)

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

  # ## Generate G data for test data (Revised)
  g_test <- scale(as.matrix(G[test_inds, 3:ncol(G)]),
                  center=mean_col_g, scale=sd_col_g)

  ## Compute P(FR | G, E)
  dup.post = testPosteriors(g_test, E_disc[test_inds], em.res)

  ## Check performance of models with N2 pairs
  RIVER.roc = roc(E_disc[test_inds_rev], dup.post$posterior[,2]) # RIVER
  GAM.roc = roc(E_disc[test_inds_rev], pp[test_inds]) # GAM

  if (verbose) {
    cat('*** AUC (GAM - genomic annotation model): ',round(GAM.roc$auc,3),'\n')
    cat('    AUC (RIVER): ',round(RIVER.roc$auc,3),'\n')
    cat('    P-value: ',format.pval(roc.test(RIVER.roc, GAM.roc)$p.value,digits=2,
                                    eps=0.001),'***\n')
    cat('\n')
  }

  list(rocRIVER=RIVER.roc, rocGAM=GAM.roc, pValue=roc.test(RIVER.roc, GAM.roc)$p.value)
}
