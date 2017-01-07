#' Application of RIVER
#'
#' \code{appRIVER} trains RIVER with all instances and computes posterior
#'         probabilities of FR for downstream analyses.
#'
#' @param G Genomic features in the data frame which contains subject IDs, gene
#'         symbols, and P number of genomic fetures of each instance per row.
#' @param E Expression outliers in the data frame which contains subject IDs,
#'         gene symbols, and outlier status of each instance per row.
#' @param pseudoc Pseudo count.
#' @param theta_init Initial values of theta.
#' @param costs Candidate penalty parameter values for L2-regularized logistic
#'         regression.
#' @param verbose Logical option for showing extra information on progress.
#'
#' @return outRIVER An object of class RIVER which contains subject IDs, gene
#'         names, posterior probabilities from GAM and RIVER, and estimated
#'         parameters from RIVER
#'
#' @section Warning: To input a vector of candidate penalty values makes
#'         \code{glmnet} faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{predict}},
#'         \code{\link{integratedEM}}, and \code{\link{testPosteriors}}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' outRIVER <- appRIVER(G, E, verbose=TRUE)
#'
#' @export

appRIVER <- function(G, E, pseudoc=50, theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  E_disc = as.numeric(E[,"Outlier"])

  g_all <- scale(as.matrix(G[,3:ncol(G)]))

  ## Search a best lambda from a multivariate logistic regression with outlier status
  ##        with 10 cross-validation
  cv.all.ll = cv.glmnet(g_all, as.vector(E_disc), lambda=costs,
                        family="binomial", alpha=0, nfolds=10) # GAM
  if (verbose) cat(' *** best lambda = ',cv.all.ll$lambda.min,' *** \n\n', sep='')

  ## Compute P(FR=1 | G)
  pp_all = predict(cv.all.ll, g_all, s="lambda.min", type="response")

  ## Train RIVER with all data for application
  em.all.res <- integratedEM(g_all, E_disc, cv.all.ll$lambda.min, cv.all.ll$glmnet.fit,
                             pseudoc, theta_init, costs, verbose)

  ## Compute P(FR | G, E)
  train.post = testPosteriors(g_all, E_disc, em.all.res)

  ## Store posterior probabilities from two models for future analyses
  list(SubjectID=G[,1], GeneName=G[,2], pFRgivenG=as.vector(pp_all),
       pFRgivenGE=train.post$posterior[,2], fitRIVER=em.all.res)
}
