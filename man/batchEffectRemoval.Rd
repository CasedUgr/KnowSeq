% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batchEffectRemoval.R
\name{batchEffectRemoval}
\alias{batchEffectRemoval}
\title{Corrects the batch effect of the data by using the selected method.}
\usage{
batchEffectRemoval(
  expressionMatrix,
  labels,
  method = "combat",
  batchGroups = c()
)
}
\arguments{
\item{expressionMatrix}{The original expression matrix to treat the batch effect.}

\item{labels}{A vector that contains the labels of the samples in expressionMatrix.}

\item{method}{The method that will be used to remove the batch effect. The possibilities are "combat" or "sva". Next release will add RUV.}

\item{batchGroups}{A numeric vector with the different known batch groups for the samples.}
}
\value{
A matrix with the batch effect corrected for combat or a model for \code{\link{DEGsExtraction}} function in the case of sva.
}
\description{
This function corrects the batch effect of the expression matrix indicated by parameter. There are two method to choose such as ComBat or SVA.
}
\examples{
dir <- system.file("extdata", package="KnowSeq")
load(paste(dir,"/expressionExample.RData",sep = ""))

batchGroups <- c(1,1,1,1,2,2,1,2,1,2)

expressionMatrixNoBatch <- batchEffectRemoval(expressionMatrix, labels, batchGroups = batchGroups)
expressionMatrixNoBatch <- batchEffectRemoval(expressionMatrix, labels, method = "sva")
}
