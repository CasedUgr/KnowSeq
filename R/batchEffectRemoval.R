#' Corrects the batch effect of the data by using the selected method.
#'
#' This function corrects the batch effect of the expression matrix indicated by parameter. There are two method to choose such as ComBat or SVA.
#'
#' @param expressionMatrix The original expression matrix to treat the batch effect.
#' @param labels A vector that contains the labels of the samples in expressionMatrix.
#' @param method The method that will be used to remove the batch effect. The possibilities are "combat" or "sva". Next release will add RUV.
#' @param batchGroups A numeric vector with the different known batch groups for the samples.
#' @return A matrix with the batch effect corrected for combat or a model for \code{\link{DEGsExtraction}} function in the case of sva.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' batchGroups <- c(1,1,1,1,2,2,1,2,1,2)
#' 
#' expressionMatrixNoBatch <- batchEffectRemoval(expressionMatrix, labels, batchGroups = batchGroups)
#' expressionMatrixNoBatch <- batchEffectRemoval(expressionMatrix, labels, method = "sva")

batchEffectRemoval <- function(expressionMatrix,labels,method = "combat", batchGroups = c()){
  
  if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
  if(typeof(expressionMatrix) != "double"){stop("The type of expression matrix elements must be double.")}
  
  if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
  
  if(!is.character(labels)){labels = as.factor(labels)}
  
  if(method == "combat"){
    
    if(length(batchGroups) == 0){stop("batchGroups parameter can not be an empty vector. The known batch groups are required to apply ComBat algorithm.")}
    
    if(!is.numeric(batchGroups)){stop("batchGroups parameter must be a numeric vector.")}
    
    if(length(batchGroups) != ncol(expressionMatrix)){stop("batchGroups parameter must be a vector with the same length than the number of columns of expressionMatrix.")}
    
    cat("Correcting batch effect by using combat method...\n")
    
    expressionMatrixBatchCorrected <- ComBat(expressionMatrix, batch = batchGroups, mod = NULL, par.prior = TRUE, mean.only = TRUE)
    
    invisible(expressionMatrixBatchCorrected)
    
  }else if(method == "sva"){
    
    cat("Calculating sva model to batch effect correction...\n")
    
    if(is.character(labels)){ labels <- as.factor(labels) }
    
    mod = model.matrix(~labels)
    mod0 = model.matrix(~1, data = labels)
    n.sv = num.sv(expressionMatrix,mod,method="leek")
    svobj = sva(expressionMatrix,mod,mod0,n.sv=n.sv)
    ndb = dim(expressionMatrix)[2]
    nmod = dim(mod)[2]
    n.sv<-svobj$n.sv
    mod1 = cbind(mod,svobj$sv)
    gammahat = (expressionMatrix %*% mod1 %*% solve(t(mod1) %*% mod1))[,(nmod+1):(nmod + svobj$n.sv)]
    expressionMatrixBatchCorrected = expressionMatrix - gammahat %*% t(svobj$sv)

    invisible(expressionMatrixBatchCorrected)
    
  }else{
    stop("The batch effect correction method selected is wrong. Please use combat or sva")
  }
  
  
}
