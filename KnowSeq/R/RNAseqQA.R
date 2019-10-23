#' RNAseqQA performs the quality analysis of an expression matrix.
#'
#' RNAseqQA performs the quality analysis of an expression matrix. This function adapts the RNA-seq data in order to allows using arrayQualityMetrics expression analysis.
#' @param expressionMatrix A matrix that contains the gene expression values.
#' @param outdir The output directory to store the report of arrayQualityMetrics
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' RNAseqQA(expressionMatrix)

RNAseqQA <- function(expressionMatrix, outdir = "RNAseqQA"){

  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}

  cat("Performing samples quality analysis...\n")

  expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]

  expressionDataset <- ExpressionSet(assayData = expressionMatrix)

  arrayQualityMetrics(expressionDataset, outdir = outdir, force = TRUE)

}
