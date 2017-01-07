#' Posterior probabilities of FR given G and E.
#'
#' \code{getFRPosteriors} computes posterior probabilities of FR (functionality
#'         of regulatory variant) given G (genomic annotations) and E (outlier
#'         status) with current estimate of beta (parameters between FR and G)
#'         and theta (parameters between FR and E).
#'
#' @param E Binary values of outlier status.
#' @param p.FR.givenG P(FR | G, beta) from \code{getFRgivenG}.
#' @param theta Current estimate of theta.
#'
#' @return P(FR | G, E, beta, theta) and probable status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' p.FR.givenG = getFRgivenG(g_all, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
#' posteriors = getFRPosteriors(E_disc, p.FR.givenG, theta=theta_init)
#'
#' @export

getFRPosteriors <- function(E, p.FR.givenG, theta) {
  p.e.FR.d1 = matrix(NA,length(E), 2)
  p.e.d1 = matrix(NA, length(E), 1)

  p.e.FR.d1 = theta[E+1,]
  p.e.d1 = rowSums(p.e.FR.d1*cbind(1.0-p.FR.givenG,p.FR.givenG))
  post = p.e.FR.d1 * c(1-p.FR.givenG, p.FR.givenG) / c(p.e.d1, p.e.d1)

  list(posterior=post, mle = max.col(post)-1)
}

#' Maximum likelihoood estimate of theta.
#'
#' \code{mleTheta} computes maximum likelihoood estimate of theta (parameters
#'         between FR (functionality of regulatory variant) and E (outlier
#'         status); Naive-Bayes).
#'
#' @param E Binary values of outlier status.
#' @param FR Soft-assignments of FR from E-step
#' @param pseudocount Pseudo count.
#'
#' @return MLE of theta
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' p.FR.givenG = getFRgivenG(g_all, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
#' posteriors = getFRPosteriors(E_disc, p.FR.givenG, theta=theta_init)
#' theta.cur = mleTheta(E_disc, posteriors$posterior, pseudoc=50)
#'
#' @export

mleTheta <- function(E, FR, pseudocount) {
  ct = matrix(NA, 2, 2)
  ct[1,1] = sum((E==0)*FR[,1])
  ct[1,2] = sum((E==0)*FR[,2])
  ct[2,1] = sum((E==1)*FR[,1])
  ct[2,2] = sum((E==1)*FR[,2])
  ct = ct + pseudocount
  return(ct/rbind(colSums(ct), colSums(ct)))
}

#' Maximum likelihoood estimate of beta.
#'
#' \code{mleBeta} computes maximum likelihoood estimate of beta (parameters
#'         between FR (functionality of regulatory variant) and G (genomic
#'         annotations); multivariate logistic regression).
#'
#' @param G Genomic features
#' @param FR Soft-assignments of FR from E-step
#' @param costs Candidate penalty parameter values for L2-regularization within
#'         logistic regression.
#'
#' @return MLE of beta
#'
#' @section Warning: To input a vector of candidate penalty values makes
#'         \code{glmnet} faster than to input a single penalty value
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' p.FR.givenG = getFRgivenG(g_all, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
#' posteriors = getFRPosteriors(E_disc, p.FR.givenG, theta=theta_init)
#' logistic.cur = mleBeta(g_all, posteriors$posterior, costs)
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @export

mleBeta <- function(G, FR, costs) {
  glmnet(G, FR, lambda=costs, family="binomial", alpha = 0)
}

#' Posterior probabilities of FR given G
#'
#' \code{getFRgivenG} computes posterior probabilities of FR (functionality of
#'         regulatory variant) given G (genomic annotations) and current estimate
#'         of beta (parameters between FR and G).
#'
#' @param G Genomic features
#' @param logistic.model Logistic regression model with current estimate of beta
#' @param lambda Selected lambda
#'
#' @return P(FR | G)
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' p.FR.givenG = getFRgivenG(g_all, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
#'
#' @seealso \code{\link{predict}}
#'
#' @export

getFRgivenG <- function(G, logistic.model, lambda) {
  predict(logistic.model, G, s=lambda, type="response")
}

#' Test posterior probabilities of FR given G and E
#'
#' \code{testPosteriors} computes posterior probabilities of FR (functionality
#'         of regulatory variant) given G (genomic annotations) and E (outlier
#'         status) with estimate of beta (parameters between FR and G) and
#'         theta (parameters between FR and E).
#'
#' @param G Genomic features
#' @param E Binary values of outlier status.
#' @param em.results Estimated parameters including beta and theta via EM and
#'         selected lambdas
#'
#' @return P(FR | G, E, beta, theta) and probable status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' em.all.res <- integratedEM(g_all, E_disc, cv.all.ll$lambda.min,
#' cv.all.ll$glmnet.fit, pseudoc=50, theta_init, costs, verbose=FALSE)
#' train.post = testPosteriors(g_all, E_disc, em.all.res)
#'
#' @seealso \code{\link{getFRgivenG}} and \code{\link{getFRPosteriors}}
#'
#' @export

testPosteriors <- function(G, E, em.results) {
  gFR = getFRgivenG(G, em.results$logistic.model, em.results$lambda)
  getFRPosteriors(E, gFR, em.results$theta)
}

#' An iterative expectation-maximization algorithm for RIVER
#'
#' \code{integratedEM} iteratively executes e-step and m-step until it
#'         converges. This is a main function of RIVER.
#'
#' @param G Genomic features.
#' @param E Binary values of outlier status.
#' @param lambda Selected lambda.
#' @param logistic.init Smart initialization of beta (parameters between
#'         FR and G) from estimate of beta with E via multivariate logistic
#'         regression.
#' @param pseudoc Pseudo count.
#' @param theta.init Initial theta (parameters between FR (functionality
#'         of regulatory variant) and E).
#' @param costs Candidate penalty parameter values for L2-regularization
#'         within logistic regression.
#' @param verbose Logical option for showing extra information on progress.
#'
#' @return Best estimate of beta and theta, final multivariate logistic
#'         regression model, and posterior probabilities of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @seealso \code{\link{getFRgivenG}}, \code{\link{getFRPosteriors}},
#'         \code{\link{mleTheta}}, \code{\link{mleBeta}},
#'         \code{\link[glmnet]{cv.glmnet}},
#'         and \url{https://github.com/ipw012/RIVERpkg}
#'
#' @examples
#' G <- simulated_features
#' E <- simulated_outliers
#' g_all <- scale(as.matrix(G[,3:ncol(G)]))
#' E_disc = as.vector(as.numeric(E[,"Outlier"]))
#' theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs,
#'         family="binomial", alpha = 0, nfolds=10)
#' em.all.res <- integratedEM(g_all, E_disc, cv.all.ll$lambda.min,
#'         cv.all.ll$glmnet.fit, pseudoc=50, theta_init, costs, verbose=FALSE)
#'
#' @export

integratedEM <- function(G, E, lambda, logistic.init, pseudoc, theta.init, costs,
                         verbose=FALSE) {
  theta.cur = theta.init
  beta.cur = logistic.init$beta[,which(logistic.init$lambda == lambda)]
  logistic.cur = logistic.init

  steps = 1
  while(TRUE) {
    if (verbose) cat(' *** RIVER: EM step ',steps,'\n',sep="")

    ## E-step:
    ## Compute expected posterior probabilities given current parameters and data
    p.FR.givenG = getFRgivenG(G, logistic.cur, lambda)
    posteriors = getFRPosteriors(E, p.FR.givenG, theta.cur)
    if (verbose) cat('     E-step: Top 10 % Threshold of expected P(FR=1 | G, E): ',
                     round(quantile(posteriors$posterior[,2], .9),4),'\n',sep='')

    ## M-step:
    ## Update theta and beta
    theta.old = theta.cur
    theta.cur = mleTheta(E, posteriors$posterior, pseudoc) # ML estimate of theta

    beta.old = beta.cur
    logistic.cur = mleBeta(G, posteriors$posterior, costs) # ML estimate of beta
    beta.cur = logistic.cur$beta[,which(logistic.cur$lambda == lambda)]

    if (verbose) cat('     M-step: norm(theta difference) = ',
                     round(norm(matrix(theta.cur)-matrix(theta.old)),4),
                     ', norm(beta difference) = ',
                     round(norm(matrix(beta.cur)-matrix(beta.old)),4),
                     " *** \n\n", sep="")

    ## Check convergence
    if ((norm(matrix(beta.cur) - matrix(beta.old)) < 1e-3) &
        (norm(theta.cur - theta.old) < 1e-3)) {
      if (verbose) cat(" ::: EM iteration is terminated (it converges within
                       the predefined tolerance, 0.001) ::: \n\n\n",sep="")
      break
    }
    steps = steps + 1
  }

  list(logistic.model=logistic.cur, beta=beta.cur, theta=theta.cur,
       posteriors=posteriors, lambda=lambda)
}
