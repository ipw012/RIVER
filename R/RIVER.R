#' RIVER: R package for an implementation of an extensible predictive model for
#'         combining whole genome sequences with molecular phenotypes to identify
#'         high-impact rare variants.
#'
#' The RIVER package provides two applicable functions of RIVER to actual
#'         datasets including genomic features and outlier status and a list of
#'         required functions for RIVER: evalRIVER, appRIVER, getFRPosteriors,
#'         mleTheta, mleBeta, getFRgivenG, testPosteriors, and integratedEM.
#'
#' @section getFuncRv_Feat: The \code{getFuncRv_Feat} computes posterior probabilities
#'         of FR (functionality of regulatory variant) given genomic features (G)
#'         and current estimate of beta (parameters between FR and G).
#' @section getFuncRvPosteriors: The \code{getFuncRvPosteriors} computes posterior
#'         probabilities of FR (functionality of regulatory variant) given genomic
#'         features (G) and outlier status (E) with current estimate of beta
#'         (parameters between FR and G) and theta (parameters between FR and
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
#' @section plotPosteriors: The \code{plotPosteriors} draw scatter plots of
#'         posterior probabilities from both RIVER GAM in terms of outliers status.
#' @section evalRIVER: The \code{evaRIVER} trains RIVER by holding out a list of
#'         individual and gene pairs having same rare variants for evaluation,
#'         computes test posterior probabilities of FR (functionality of regulatory
#'         variant) for 1st individual, and compares them with outlier status of
#'         2nd individual from the list.
#' @section appRIVER: The \code{appRIVER} trains RIVER with all instances and
#'         computes posterior probabilities of FR (functionality of regulatory
#'         variant) for downstream analyses.
#'
#' @source \url{https://github.com/ipw012/RIVER}
#'
#' @docType package
#' @name RIVER
#' @import dplyr pROC ggplot2
#' @rawNamespace import(glmnet, except = "auc")
#' @importFrom data.table fread
#' @importFrom Biobase ExpressionSet exprs sampleNames
#' @importFrom methods new
#' @importFrom graphics abline legend lines par plot title
#' @importFrom stats predict quantile sd
#' @importFrom glmnet cv.glmnet
#' @importFrom utils read.table
NULL
