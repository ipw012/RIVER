#' Application of RIVER
#'
#' \code{appRIVER} trains RIVER with all instances and computes posterior
#'         probabilities of FR for downstream analyses.
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
#' @return A list which contains subject IDs, gene names, posterior
#'         probabilities from GAM and RIVER, and estimated parameters from
#'         RIVER with used hyperparameters.
#'
#' @section Warning: To input a vector of candidate penalty values makes
#'         \code{glmnet} faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{predict}},
#'         \code{\link{integratedEM}}, \code{\link{testPosteriors}},
#'         \code{\link{getData}}, \code{\link[Biobase]{exprs}}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVERpkg"), ZscoreThrd=1.5)
#' postprobs <- appRIVER(dataInput, verbose=TRUE)
#'
#' @export

appRIVER <- function(dataInput, pseudoc=50, theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Extract required data for evaRIVER
  GAll = t(exprs(dataInput)) # all genomic features (G)
  EAll = as.numeric(unlist(dataInput$Outlier))-1 # all outlier status (E)

  ## Search a best lambda from a multivariate logistic regression with outlier status
  ##        with 10 cross-validation
  cv.all.ll = cv.glmnet(GAll, as.vector(EAll), lambda=costs,
                        family="binomial", alpha=0, nfolds=10) # GAM
  if (verbose) cat(' *** best lambda = ',cv.all.ll$lambda.min,' *** \n\n', sep='')

  ## Compute P(FR=1 | G)
  pp_all = predict(cv.all.ll, GAll, s="lambda.min", type="response")

  ## Train RIVER with all data for application
  em.all.res <- integratedEM(GAll, EAll, cv.all.ll$lambda.min, cv.all.ll$glmnet.fit,
                             pseudoc, theta_init, costs, verbose)

  ## Compute P(FR | G, E)
  train.post = testPosteriors(GAll, EAll, em.all.res)

  ## Output: postprobs
  postprobs=list(indiv_name=unlist(strsplit(rownames(GAll),":"))[seq(1,2*nrow(GAll),by=2)],
       gene_name=unlist(strsplit(rownames(GAll),":"))[seq(2,2*nrow(GAll),by=2)],
       RIVER_posterior=train.post$posterior[,2],
       GAM_posterior=as.numeric(pp_all),
       fitRIVER=em.all.res)
  class(postprobs) = "appl"
  return(postprobs)
}
