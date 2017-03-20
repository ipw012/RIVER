#' Application of RIVER
#'
#' \code{appRIVER} trains RIVER with all instances and computes posterior
#'         probabilities of FR for downstream analyses.
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
#'         package = "RIVER"), ZscoreThrd=1.5)
#' postprobs <- appRIVER(dataInput, verbose=TRUE)
#'
#' @export

appRIVER <- function(dataInput, pseudoc=50, theta_init=matrix(c(.99, .01, .3, .7), nrow=2),
                     costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4), verbose=FALSE) {
  ## Extract required data for evaRIVER
  FeatAll = t(exprs(dataInput)) # all genomic features (G)
  OutAll = as.numeric(unlist(dataInput$Outlier))-1 # all outlier status (E)

  ## Search a best lambda from a multivariate logistic regression with outlier status
  ##        with 10 cross-validation
  logisticAllCV = cv.glmnet(FeatAll, as.vector(OutAll), lambda=costs,
                        family="binomial", alpha=0, nfolds=10) # GAM
  if (verbose) cat(' *** best lambda = ',logisticAllCV$lambda.min,' *** \n\n', sep='')

  ## Compute P(FR=1 | G)
  postporbGAM = predict(logisticAllCV, FeatAll, s="lambda.min", type="response")

  ## Train RIVER with all data for application
  emModelAll <- integratedEM(FeatAll, OutAll, logisticAllCV$lambda.min, logisticAllCV$glmnet.fit,
                             pseudoc, theta_init, costs, verbose)

  ## Compute P(FR | G, E)
  postprobRIVER = testPosteriors(FeatAll, OutAll, emModelAll)

  ## Output: postprobs
  postprobs=list(indiv_name=unlist(strsplit(rownames(FeatAll),":"))[seq(1,2*nrow(FeatAll),by=2)],
       gene_name=unlist(strsplit(rownames(FeatAll),":"))[seq(2,2*nrow(FeatAll),by=2)],
       RIVER_posterior=postprobRIVER$posterior[,2],
       GAM_posterior=as.numeric(postporbGAM),
       fitRIVER=emModelAll)
  class(postprobs) = "appl"
  return(postprobs)
}
