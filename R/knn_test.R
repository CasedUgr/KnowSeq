#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset.
#'
#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset. An optimization of the k neighbours is done at the start of the process.
#'
#' @param test The test parameter is an expression matrix or data.frame that contains the test dataset with the genes in the columns and the samples in the rows.
#' @param labelsTest A vector or factor that contains the test labels for each of the samples in the test object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param bestK Best K selected during the training phase.
#' @return A list that contains four objects. The confusion matrix, the accuracy, the sensitibity and the specificity for each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' testMatrix <- t(DEGsMatrix)[c(5,10),]
#' testLabels <- labels[c(5,10)]
#' bestK <- 3 # the one that has been selectede
#' results_test_knn <- knn_test(testMatrix, testLabels, rownames(DEGsMatrix)[1:10], bestK)

knn_test <-function(test,labelsTest,vars_selected, bestK){

  if(!is.character(labelsTest)  && !is.factor(labelsTest)){stop("The class of the labelsTest parameter must be character vector or factor.")}
  if(is.character(labelsTest)){ labelsTest <- as.factor(labelsTest) }

  if(!is.data.frame(test) && !is.matrix(test)){

    stop("The test argument must be a dataframe or a matrix.")

  }

  if(dim(test)[1] != length(labelsTest)){

    stop("The length of the rows of the argument test must be the same than the length of the lablesTest. Please, ensures that the rows are the samples and the columns are the variables.")

  }

  test <- as.data.frame(apply(test,2,as.double))
  test <- test[,vars_selected]
  
  test = vapply(test, function(x){ 
    max = max(x)
    min = min(x)
    x = ((x-min)/(max-min))*2-1}, double(nrow(test)))
  
  test <- as.data.frame(test)

  accVector <- double()
  sensVector <- double()
  specVector <- double()
  cfMatList  <- list()
  
  # Firstly with one variable
  cat(paste("Testing with ", 1," variables...\n",sep=""))
  predicts<-knn(train[, 1, drop=FALSE],test[, 1, drop=FALSE],labelsTrain,k=bestK)
  
  cfMat<-confusionMatrix(predicts,labelsTest)
  acc<-confusionMatrix(predicts,labelsTest)$overall[[1]]
  sens<-confusionMatrix(predicts,labelsTest)$byClass[[1]]
  spec<-confusionMatrix(predicts,labelsTest)$byClass[[2]]
  
  cfMatList[[1]] <- cfMat
  accVector[1] <- acc
  sensVector[1] <- sens
  specVector[1] <- spec
  
  for(i in c(2:dim(train)[2])){

    cat(paste("Testing with ", i," variables...\n",sep=""))
    predicts<-knn(train[,seq_len(i)],test[,seq_len(i)],labelsTrain,k=bestK)

    cfMat<-confusionMatrix(predicts,labelsTest)
    acc<-confusionMatrix(predicts,labelsTest)$overall[[1]]
    sens<-confusionMatrix(predicts,labelsTest)$byClass[[1]]
    spec<-confusionMatrix(predicts,labelsTest)$byClass[[2]]

    cfMatList[[i]] <- cfMat
    accVector[i] <- acc
    sensVector[i] <- sens
    specVector[i] <- spec

  }

  cat("Classification done successfully!\n")
  names(accVector) <- vars_selected
  names(sensVector) <- vars_selected
  names(specVector) <- vars_selected

  results <- list(cfMatList,accVector,sensVector,specVector)
  names(results) <- c("cfMats","accVector","sensVector","specVector")
  invisible(results)

}
