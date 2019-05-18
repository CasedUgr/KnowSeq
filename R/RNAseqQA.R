#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function adapts the RNA-seq data in order to allows using arrayQualityMetrics expression analysis.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @return Nothing to return.
#' @examples
#' downloadPublicSeries(c("GSE74251","GSE81593"))
#'
#' GSE74251 <- read.csv("ReferenceFiles/GSE74251.csv")
#' GSE81593 <- read.csv("ReferenceFiles/GSE81593.csv")
#'
#' GSE74251 <- GSE74251[1:5,]
#' GSE81593 <- GSE81593[8:12,]
#'
#' dir <- system.file("extdata", package="KnowSeq")
#'
#' Run <- GSE74251$Run
#' Path <- paste(dir,"/countFiles/",GSE74251$Run,sep = "")
#' Class <- rep("Tumor", length(GSE74251$Run))
#' GSE74251CountsInfo <-  data.frame(Run = Run, Path = Path, Class = Class)
#'
#' Run <- GSE81593$Run
#' Path <- paste(dir,"/countFiles/",GSE81593$Run,sep = "")
#' Class <- rep("Control", length(GSE81593$Run))
#' GSE81593CountsInfo <-  data.frame(Run = Run, Path = Path, Class = Class)
#'
#' mergedCountsInfo <- rbind(GSE74251CountsInfo, GSE81593CountsInfo)
#'
#' write.csv(mergedCountsInfo, file = "ReferenceFiles/mergedCountsInfo.csv")
#'
#' countsInformation <- countsToMatrix("ReferenceFiles/mergedCountsInfo.csv")
#'
#' countsMatrix <- countsInformation$countsMatrix
#' labels <- countsInformation$labels
#'
#' myAnnotation <- getAnnotationFromEnsembl(rownames(countsMatrix),referenceGenome=37)
#'
#' expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation, genesNames = TRUE)
#'
#' RNAseqQA(expressionMatrix)

RNAseqQA <- function(expressionMatrix, outdir = "RNAseqQA"){

  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}

  cat("Performing samples quality analysis...\n")

  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]

  expressionDataset <- ExpressionSet(assayData = expressionMatrix)

  arrayQualityMetrics(expressionDataset, outdir = outdir, force = TRUE)

}
