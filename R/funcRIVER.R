#' Posterior probabilities of FR given G and E.
#'
#' \code{getFuncRvPosteriors} computes posterior probabilities of functionality
#'         of regulatory variant (FR) given genomic features (G) and outlier
#'         status (E) with current estimate of beta (parameters between FR and G)
#'         and theta (parameters between FR and E).
#'
#' @param Out Binary values of outlier status (E).
#' @param probFuncRvFeat probabilities of FR given genomic features and estimated
#'         beta, P(FR | G, beta), from \code{getFuncRvFeat}.
#' @param theta Current estimate of theta.
#'
#' @return posterior probabilities of FR (P(FR | G, E, beta, theta)) and probable
#'         status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init<-matrix(c(.99, .01, .3, .7), nrow=2)
#' costs<-c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logisticAllCV$glmnet.fit,
#'         logisticAllCV$lambda.min)
#' posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta=theta.init)
#'
#' @export

getFuncRvPosteriors <- function(Out, probFuncRvFeat, theta) {
  probOut_FuncRv <- matrix(NA,length(Out), 2)
  probOut <- matrix(NA, length(Out), 1)

  probOut_FuncRv <- theta[Out+1,]
  probOut <-
    rowSums(probOut_FuncRv*cbind(1.0-probFuncRvFeat,probFuncRvFeat))
  post <-
    probOut_FuncRv*c(1-probFuncRvFeat,probFuncRvFeat)/c(probOut, probOut)

  list(posterior=post, mle = max.col(post)-1)
}

#' Maximum likelihoood estimate of theta.
#'
#' \code{mleTheta} computes maximum likelihoood estimate of theta (parameters
#'         between FR (functionality of regulatory variant) and E (outlier
#'         status); Naive-Bayes).
#'
#' @param Out Binary values of outlier status (E).
#' @param FuncRv Soft-assignments of FR from E-step
#' @param pseudocount Pseudo count.
#'
#' @return MLE of theta
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init <- matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logisticAllCV$glmnet.fit, logisticAllCV$lambda.min)
#' posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta=theta.init)
#' thetaCur <- mleTheta(Out, FuncRv=posteriors$posterior, pseudoc=50)
#'
#' @export

mleTheta <- function(Out, FuncRv, pseudocount) {
  countOut <- matrix(NA, 2, 2)
  countOut[1,1] <- sum((Out==0)*FuncRv[,1])
  countOut[1,2] <- sum((Out==0)*FuncRv[,2])
  countOut[2,1] <- sum((Out==1)*FuncRv[,1])
  countOut[2,2] <- sum((Out==1)*FuncRv[,2])
  countOut <- countOut + pseudocount
  return(countOut/rbind(colSums(countOut), colSums(countOut)))
}

#' Maximum likelihoood estimate of beta.
#'
#' \code{mleBeta} computes maximum likelihoood estimate of beta (parameters
#'         between FR (functionality of regulatory variant) and G (genomic
#'         annotations); multivariate logistic regression).
#'
#' @param Feat Genomic features (G)
#' @param FuncRv Soft-assignments of FR from E-step
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
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init <- matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha=0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logisticAllCV$glmnet.fit, logisticAllCV$lambda.min)
#' posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta=theta.init)
#' logistic.cur <- mleBeta(Feat, FuncRv=posteriors$posterior, costs)
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @export

mleBeta <- function(Feat, FuncRv, costs) {
  glmnet(Feat, FuncRv, lambda=costs, family="binomial", alpha = 0)
}

#' Posterior probabilities of FR given G
#'
#' \code{getFuncRvFeat} computes posterior probabilities of FR (functionality of
#'         regulatory variant) given G (genomic features) and current estimate
#'         of beta (parameters between FR and G).
#'
#' @param Feat Genomic features (G)
#' @param logistic.model Logistic regression model with current estimate of beta
#' @param lambda Selected lambda
#'
#' @return probabilities of FR given genomic features, P(FR | G)
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' probFuncRvFeat <- getFuncRvFeat(Feat, logistic.model=logisticAllCV$glmnet.fit,
#'         lambda=logisticAllCV$lambda.min)
#'
#' @seealso \code{\link{predict}}
#'
#' @export

getFuncRvFeat <- function(Feat, logistic.model, lambda) {
  predict(logistic.model, Feat, s=lambda, type="response")
}

#' Test posterior probabilities of FR given G and E
#'
#' \code{testPosteriors} computes posterior probabilities of FR (functionality
#'         of regulatory variant) given G (genomic annotations) and E (outlier
#'         status) with estimate of beta (parameters between FR and G) and
#'         theta (parameters between FR and E).
#'
#' @param Feat Genomic features (G)
#' @param Out Binary values of outlier status (E).
#' @param emModel Estimated parameters including beta and theta via EM and
#'         selected lambdas
#'
#' @return test posterior probabilities of FR given new outlier status (E)
#'         and genomic features (G), P(FR | G, E, beta, theta), and probable
#'         status of FR.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init <- matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' emModelAll <- integratedEM(Feat, Out, logisticAllCV$lambda.min, logisticAllCV$glmnet.fit,
#'         pseudoc=50, theta.init, costs, verbose=FALSE)
#' trainedpost <- testPosteriors(Feat, Out, emModel=emModelAll)
#'
#' @seealso \code{\link{getFuncRvFeat}} and \code{\link{getFuncRvPosteriors}}
#'
#' @export

testPosteriors <- function(Feat, Out, emModel) {
  probFuncRvFeat <- getFuncRvFeat(Feat, emModel$logistic.model, emModel$lambda)
  getFuncRvPosteriors(Out, probFuncRvFeat, emModel$theta)
}

#' An iterative expectation-maximization algorithm for RIVER
#'
#' \code{integratedEM} iteratively executes e-step and m-step until it
#'         converges. This is a main function of RIVER.
#'
#' @param Feat Genomic features (G).
#' @param Out Binary values of outlier status (E).
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
#' @seealso \code{\link{getFuncRvFeat}}, \code{\link{getFuncRvPosteriors}},
#'         \code{\link{mleTheta}}, \code{\link{mleBeta}},
#'         \code{\link[glmnet]{cv.glmnet}},
#'         and \url{https://github.com/ipw012/RIVER}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' Feat <- scale(t(Biobase::exprs(dataInput))) # genomic features (G)
#' Out <- as.vector(as.numeric(unlist(dataInput$Outlier))-1) # outlier status (E)
#' theta.init=matrix(c(.99, .01, .3, .7), nrow=2)
#' costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)
#' logisticAllCV <- glmnet::cv.glmnet(Feat, Out, lambda=costs, family="binomial",
#'         alpha = 0, nfolds=10)
#' emModelAll <- integratedEM(Feat, Out, lambda=logisticAllCV$lambda.min,
#'         logistic.init=logisticAllCV$glmnet.fit, pseudoc=50, theta=theta.init,
#'         costs, verbose=FALSE)
#'
#' @export

integratedEM <- function(Feat, Out, lambda, logistic.init,
                         pseudoc, theta.init, costs,
                         verbose=FALSE){
  theta.cur <- theta.init
  beta.cur <- logistic.init$beta[,which(logistic.init$lambda == lambda)]
  logistic.cur <- logistic.init

  steps <- 1
  maxIter <- 1000
  converged <- 0
  for (iter in 1:maxIter) {
    if (verbose) {
      cat(' *** RIVER: EM step ',steps,'\n',sep="")
    }

    ## E-step:
    ## Compute expected posterior probabilities
    ##           given current parameters and data
    probFuncRvFeat <- getFuncRvFeat(Feat, logistic.cur, lambda)
    posteriors <- getFuncRvPosteriors(Out, probFuncRvFeat, theta.cur)
    if (verbose) {
      cat('     E-step: Top 10 % Threshold of expected P(FR=1 | G, E): ',
          round(quantile(posteriors$posterior[,2], .9),4),'\n',sep='')
    }

    ## M-step:
    ## Update theta and beta
    theta.old <- theta.cur
    # ML estimate of theta
    theta.cur <- mleTheta(Out, posteriors$posterior, pseudoc)

    beta.old <- beta.cur
    # ML estimate of beta
    logistic.cur <- mleBeta(Feat, posteriors$posterior, costs)
    beta.cur <- logistic.cur$beta[,which(logistic.cur$lambda == lambda)]

    if (verbose) {
      cat('     M-step: norm(theta difference) = ',
          round(norm(matrix(theta.cur)-matrix(theta.old)),4),
          ', norm(beta difference) = ',
          round(norm(matrix(beta.cur)-matrix(beta.old)),4),
          " *** \n\n", sep="")
    }

    ## Check convergence
    if ((norm(matrix(beta.cur) - matrix(beta.old)) < 1e-3) &
        (norm(theta.cur - theta.old) < 1e-3)) {
      converged <- 1
      break
    }
    steps <- steps + 1
  }

  if (converged == 1) {
    cat(" ::: EM iteration is terminated since it converges within a
        predefined tolerance (0.001) ::: \n\n\n",sep="")
  } else if ((converged == 0) && (iter == maxIter)) {
    cat(" ::: EM iteration is terminated since it reaches a
        predefined maximum value (1000) ::: \n\n\n",sep="")
  }

  list(logistic.model=logistic.cur, beta=beta.cur, theta=theta.cur,
       posteriors=posteriors, lambda=lambda)
}
