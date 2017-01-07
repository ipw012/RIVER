#' RIVERpkg: R package for an implementation of an extensible predictive model for
#'         combining whole genome sequences with molecular phenotypes to identify
#'         high-impact rare variants.
#'
#' The RIVER package provides two applicable functions of RIVER to actual
#'         datasets including genomic features and outlier status and a list of
#'         required functions for RIVER: evalRIVER, appRIVER, getFRPosteriors,
#'         mleTheta, mleBeta, getFRgivenG, testPosteriors, and integratedEM.
#'
#' @section getFRgivenG: The \code{getFRgivenG} computes posterior probabilities
#'         of FR (functionality of regulatory variant) given G (genomic
#'         annotations) and current estimate of beta (parameters between FR and G).
#' @section getFRPosteriors: The \code{getFRPosteriors} computes posterior
#'         probabilities of FR (functionality of regulatory variant) given G
#'         (genomic annotations) and E (outlier status) with current estimate of
#'         beta (parameters between FR and G) and theta (parameters between FR and
#'         E).
#' @section mleBeta: The \code{mleBeta} computes maximum likelihoood estimate of
#'         beta (parameters between FR and G; multivariate logistic regression).
#' @section mleTheta: The \code{mleTheta} computes maximum likelihoood estimate
#'         of theta (parameters between FR and E; Naive-Bayes).
#' @section testPosteriors: The \code{testPosteriors} computes posterior
#'         probabilities of FR (functionality of regulatory variant) given G
#'         (genomic annotations) and E (outlier status) with estimate of beta
#'         (parameters between FR and G) and theta (parameters between FR and E).
#' @section integratedEM: The \code{integratedEM} iteratively executes e-step
#'         and m-step until it converges. This is a main function of RIVER.
#' @section evalRIVER: The \code{evaRIVER} trains RIVER by holding out a list of
#'         individual and gene pairs having same rare variants for evaluation,
#'         computes test posterior probabilities of FR (functionality of regulatory
#'         variant) for 1st individual, and compares them with outlier status of
#'         2nd individual from the list.
#' @section appRIVER: The \code{appRIVER} trains RIVER with all instances and
#'         computes posterior probabilities of FR (functionality of regulatory
#'         variant) for downstream analyses.
#'
#' @source \url{https://github.com/ipw012/RIVERpkg}
#'
#' @docType package
#' @name RIVERpkg
#' @import dplyr pROC ggplot2
#' @rawNamespace import(glmnet, except = "auc")
#' @importFrom graphics abline legend lines par plot title
#' @importFrom stats predict quantile sd
#' @importFrom glmnet cv.glmnet
NULL
