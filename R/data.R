#' Continuous values of 18 simulated genomic features.
#'
#' A simulated dataset having scaled values of 18 distinct simulated genomic
#'         features for 6122 instances. Note that simulated_features and
#'         simulated_outliers are matched in the same order.
#'
#' @format A data frame with 6122 rows and 20 columns including subject ID,
#'         gene name, and 18 features:
#'   \describe{
#'   \item{SubjectID}{(string) - subject IDs of 115 distinct individuals}
#'   \item{GeneName}{(string) - 1631 distinct gene names}
#'   \item{Feature1}{scaled values of genomic features 1}
#'   \item{Feature2}{scaled values of genomic features 2}
#'   \item{Feature3}{scaled values of genomic features 3}
#'   \item{Feature4}{scaled values of genomic features 4}
#'   \item{Feature5}{scaled values of genomic features 5}
#'   \item{Feature6}{scaled values of genomic features 6}
#'   \item{Feature7}{scaled values of genomic features 7}
#'   \item{Feature8}{scaled values of genomic features 8}
#'   \item{Feature9}{scaled values of genomic features 9}
#'   \item{Feature10}{scaled values of genomic features 10}
#'   \item{Feature11}{scaled values of genomic features 11}
#'   \item{Feature12}{scaled values of genomic features 12}
#'   \item{Feature13}{scaled values of genomic features 13}
#'   \item{Feature14}{scaled values of genomic features 14}
#'   \item{Feature15}{scaled values of genomic features 15}
#'   \item{Feature16}{scaled values of genomic features 16}
#'   \item{Feature17}{scaled values of genomic features 17}
#'   \item{Feature18}{scaled values of genomic features 18}}
#' @return A simulated feature data in a data.frame format (6122 x 20)
#'
#' @source \url{https://github.com/ipw012/RIVERpkg}
"simulated_features"

#' Binary values of gene expression outlier status.
#'
#' A simulated dataset  having binary values of gene expression outlier status
#'         for 6122 instances. Note that simulated_features and simulated_outliers
#'         are matched in the same order.
#'
#' @format A data frame with 6122 rows and 3 columns including subject ID, gene
#'         name, and outlier status:
#'   \describe{
#'   \item{SubjectID}{subject IDs of 115 distinct individuals}
#'   \item{GeneName}{1631 distinct gene names}
#'   \item{Outlier}{binary values of outlier status}}
#' @return A simulated outlier data in a data.frame format (6122 x 3)
#'
#' @source \url{https://github.com/ipw012/RIVERpkg}
"simulated_outliers"
