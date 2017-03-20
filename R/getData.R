#' Get coordinated data from a compressed file for RIVER
#'
#' \code{getData} extracts genomic features, z-scores of gene expression, and
#'         N2 pairs having same rare variants from an imported compressed data,
#'         computes outlier status from z-scores given a z-score threshold and
#'         coordinates the genomic features, outlier status, and a list of N2
#'         pairs into ExpressionSet class having standardized data structure.
#'
#' @param filename A full path of a compressed input file that consists of all
#'         samples in rows and subject ID, gene name, genomic features, z-scores
#'         of corresponding gene expression, and a list of N2 pairs in columns
#'         from left to right. In N2 pairs, samples not paired with othersamples
#'         have NA while two samples sharing same rare variant near a gene have
#'         same pre-assigend integers.
#' @param ZscoreThrd A |Z-score| threshold for defining outlier status of samples
#'
#' @return dataInput An object of ExpressionSet class which contains input data
#'         required for all functions in RIVER including genomic features,
#'         outlier status, and N2 pairs.
#'
#' @author Yungil Kim, \email{ipw012@@gmail.com}
#' @seealso \code{\link[data.table]{fread}}, \code{\link[Biobase]{ExpressionSet}},
#'
#' @examples
#' InputData <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
#'         package = "RIVER"), ZscoreThrd=1.5)
#'
#' @export

getData <- function(filename=system.file("extdata", "simulation_RIVER.gz",
                                         package = "RIVER"), ZscoreThrd=1.5) {
  # Only for linux users
  # expData = as.data.frame(fread(paste("zcat ", filename, sep=""),
  #                               sep='\t', header = TRUE, na.strings = "NA"))
  expData = read.table(gzfile(filename), header = TRUE)

  Feat = expData[,3:(ncol(expData)-2)] # genomic features
  rownames(Feat) = paste(expData[,"SubjectID"], ":",
                      expData[,"GeneName"],sep="") # sample name as SubjectID:GeneName
  Feat = as.matrix(t(Feat)) # feature x sample

  # outlier status, N2 pairs
  pData = data.frame(Outlier=factor(ifelse(abs(expData[,"Zscore"])>=ZscoreThrd,1,0),
                                    levels=c(0,1)),
                     N2pair=factor(expData[,"N2pair"],
                                   levels=unique(expData[,"N2pair"])))
  rownames(pData) = paste(expData[,"SubjectID"],":",expData[,"GeneName"],sep="")

  # descrition of outlier status and N2 pairs
  metadata = data.frame(labelDescription=c("Outlier status based on Z-scores",
                                           "Pairs of samples having same rare variants"),
                        row.names=c("Outlier","N2pair"))
  phenoData = new("AnnotatedDataFrame",
                  data=pData, varMetadata=metadata)
  dataInput = ExpressionSet(assayData=Feat, phenoData=phenoData)
  return(dataInput)
}
