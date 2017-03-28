#' Draw scatter plots of posterior probabilities from both RIVER GAM in terms
#'         of outlier status.
#'
#' \code{plotPosteriors} draws scatter plots of posterior probabilities from
#'         both RIVER GAM (genomic annotation model) in terms of outlier
#'         status.
#'
#' @param postprobs Output of \code{evaRIVER}, which provides test posterior
#'         probabilities from both RIVER and GAM for all instances.
#' @param outliers Outlier status of examples
#'
#' @return A figure of posteriors from RIVER (y-axis) and GAM (x-axis) models
#'         for ouliters and non-outliers separately
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#' postprobs <- appRIVER(dataInput)
#' plotPosteriors(postprobs, outliers=as.numeric(unlist(dataInput$Outlier))-1)
#'
#' @export

plotPosteriors <- function(postprobs, outliers) {
  probFuncRv_Feat <- probFuncRv_FeatOut <- Outliers <- NULL

  par(mar=c(6.1, 6.1, 4.1, 4.1))
  dat <- data.frame(probFuncRv_Feat=postprobs$GAM_posterior,
                   probFuncRv_FeatOut=postprobs$RIVER_posterior,
                   outliers=as.factor(outliers))

  ggplot(dat, aes_string(x="probFuncRv_Feat",
                         y="probFuncRv_FeatOut",
                         color="outliers")) +
    geom_point(shape=1, size=4) +
    geom_abline(intercept=0, slope=1, color="darkgray",
                size=1, linetype=2) +
    theme_bw() + xlab("P( FR | G)") + ylab("P( FR | G, E)") +
    scale_color_manual(values=c("dodgerblue","mediumpurple"),
                       name="", breaks=c("0", "1"),
                       labels=c("Non-outlier", "Outlier")) +
    theme(axis.title = element_text(size=18),
          axis.text  = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.15,0.9),
          legend.text = element_text(size=16),
          legend.key = element_blank())
}
