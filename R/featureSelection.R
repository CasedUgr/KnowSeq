#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest.
#'
#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest. Furthermore, the ranking is returned and can be used as input of the parameter vars_selected in the machine learning functions.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each samples in data parameter.
#' @param vars_selected The genes selected to use in the feature selection process. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes.
#' @param mode The algorithm used to calculate the genes ranking. The possibilities are two: mrmr and rf.
#' @return A vector that contains the ranking of genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' featureRanking <- featureSelection(t(DEGsMatrix),labels,rownames(DEGsMatrix))

featureSelection <-function(data,labels,vars_selected,mode="mrmr"){

  if(!is.data.frame(data) && !is.matrix(data)){

    stop("The data argument must be a dataframe or a matrix.")

  }
  if(dim(data)[1] != length(labels)){

    stop("The length of the rows of the argument data must be the same than the length of the labels. Please, ensures that the rows are the samples and the columns are the variables.")

  }

  if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
  if(is.character(labels)){ labels <- as.factor(labels) }

  if(length(vars_selected)[1] == 0 || is.null(vars_selected)){

    stop("The vars_selected parameter is empty! Please, provide a right vars list.")

  }

  data <- as.data.frame(apply(data,2,as.double))

  if(mode == "mrmr"){

    cat("Calculating the ranking of the most relevant genes by using mRMR algorithm...\n")

    mrmrRanking <- MRMR(data[,vars_selected],labels, k=ncol(data[,vars_selected]))

    cat("mRMR ranking: ")
    cat(vars_selected[mrmrRanking$selection])
    cat("\n")

    return(mrmrRanking$selection)

  }else if(mode == "rf"){

    cat("Calculating the ranking of the most relevant genes by using Random Forest algorithm...\n")

    rfRanking <- randomForest(data[,vars_selected], labels, importance=TRUE,proximity=TRUE)
    rfRanking <- rfRanking$importance[order(rfRanking$importance[,3],decreasing = TRUE),]
    
    cat("Random Forest ranking: ")
    cat(rownames(rfRanking))
    cat("\n")

    return(rownames(rfRanking))

  }else{
    stop("The mode is unrecognized, please use mrmr or rf.")
  }
}
