#' Corrects the batch effect of the data by using the selected method.
#'
#' This function corrects the batch effect of the expression matrix indicated by parameter. There are two method to choose such as ComBat or SVA.
#'
#' @param expressionMatrix The original expression matrix to treat the batch effect.
#' @param labels A vector that contains the labels of the samples in expressionMatrix.
#' @param method The method that will be used to remove the batch effect. The possibilities are "combat" or "sva". Next release will add RUV.
#' @param clusters The number of clusters intrinsic to the expression matrix data which could means different batches. The optimal number of clusters in the expression matrix can be calculated by calling the function \code{\link{dataPlot}}, with the parameter mode equal to "optimalClusters". This parameter is only required when the user selects the combat method.
#' @return A matrix with the batch effect corrected for combat or a model for \code{\link{limmaDEGsExtraction}} function in the case of sva.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' expressionMatrixNoBatch <- batchEffectRemoval(expressionMatrix, labels, clusters = 4)
#' svaMod <- batchEffectRemoval(expressionMatrix, labels, method = "sva")

batchEffectRemoval <- function(expressionMatrix,labels,method = "combat", clusters = 2){
  
  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  if(typeof(expressionMatrix) != "double"){stop("The type of expression matrix elements must be double.")}

  if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}

  if(!is.numeric(clusters)){stop("The class of clusters parameter must be numeric.")}

  if(clusters < 2){stop("Clusters parameter must be greater or equal than 2.")}

  if(!is.character(labels)){labels = as.factor(labels)}

  if(method == "combat"){

   cat("Correcting batch effect by using combat method...\n")

   km.res <- kmeans(t(expressionMatrix), clusters, nstart = 25)
   expressionMatrixBatchCorrected <- ComBat(expressionMatrix,batch = km.res$cluster, mod = NULL, par.prior = TRUE, mean.only = TRUE)



   return(expressionMatrixBatchCorrected)

  }else if(method == "sva"){

    cat("Calculating sva model to batch effect correction...\n")

    if(is.character(labels)){ labels <- as.factor(labels) }

    mod = model.matrix(~labels)
    mod0 = model.matrix(~1, data = labels)
    n.sv = num.sv(expressionMatrix,mod,method="leek")
    svobj = sva(expressionMatrix,mod,mod0,n.sv=n.sv)
    pValues = f.pvalue(expressionMatrix,mod,mod0)
    qValues = p.adjust(pValues,method="BH")
    modSv = cbind(mod,svobj$sv)
    mod0Sv = cbind(mod0,svobj$sv)
    pValuesSv = f.pvalue(expressionMatrix,modSv,mod0Sv)
    qValuesSv = p.adjust(pValuesSv,method="BH")

    return(modSv)

  }else{
    stop("The batch effect correction method selected is wrong. Please use combat or sva")
  }


}
