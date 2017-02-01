#' Draw AUC curves from both RIVER and GAM.
#'
#' \code{plotAUC} plots AUC curves from both RIVER and GAM (genomic annotation
#'         model).
#'
#' @param evaROC Output of \code{evaRIVER}, an S4 object of class evaRIVER which
#'         include two AUC values from RIVER and GAM, computed specificities and
#'         sensitivities from two models, and P-value of comparing the two AUC
#'         values.
#'
#' @return AUC figure
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#'
#' @examples
#' data <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVERpkg"), ZscoreThrd=1.5)
#' evaROC <- evaRIVER(data, verbose=TRUE)
#' plotAUC(evaROC)
#'
#' @export

plotAUC <-function(evaROC){
  par(mar=c(6.1, 6.1, 4.1, 4.1))
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="False positive rate",
       ylab="True positive rate", cex.axis=1.3, cex.lab=1.6)
  abline(0, 1, col="gray")
  lines(1-evaROC$RIVER_spec, evaROC$RIVER_sens, type="s", col='dodgerblue', lwd=2)
  lines(1-evaROC$GAM_spec, evaROC$GAM_sens, type="s", col='mediumpurple', lwd=2)
  legend(0.7,0.2,c("RIVER","GAM"), lty=c(1,1), lwd=c(2,2),
         col=c("dodgerblue","mediumpurple"), cex=1.2, pt.cex=1.2, bty="n")
  title(main=paste("AUC: RIVER = ", round(evaROC$RIVER_auc,3), ", GAM = ",
                   round(evaROC$GAM_auc,3), ", P = ",
                   format.pval(evaROC$pvalue, digits=2, eps=0.001),sep=""))
}

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
#'         package = "RIVERpkg"), ZscoreThrd=1.5)
#' postprobs <- appRIVER(dataInput)
#' plotPosteriors(postprobs, outliers=as.numeric(unlist(dataInput$Outlier))-1)
#'
#' @export

plotPosteriors <- function(postprobs, outliers) {
  pFRgivenG <- pFRgivenGE <- Outliers <- NULL

  par(mar=c(6.1, 6.1, 4.1, 4.1))
  dat = data.frame(pFRgivenG=postprobs$GAM_posterior,
                   pFRgivenGE=postprobs$RIVER_posterior,
                   outliers=as.factor(outliers))
  colnames(dat) = c("pFRgivenG","pFRgivenGE","Outliers")

  ggplot(dat, aes(x=pFRgivenG, y=pFRgivenGE, color=Outliers)) +
    geom_point(shape=1, size=4) +
    geom_abline(intercept=0,slope=1,color="darkgray",size=1,linetype=2) +
    theme_bw() + xlab("P( FR | G)") + ylab("P( FR | G, E)") +
    scale_color_manual(values=c("dodgerblue","mediumpurple"),
                       name="",
                       breaks=c("0", "1"),
                       labels=c("Non-outlier", "Outlier")) +
    theme(axis.title.x = element_text(size=18),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=18),
          axis.text.y  = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.15,0.9),
          legend.text = element_text(size=16),
          legend.key = element_blank())
}
