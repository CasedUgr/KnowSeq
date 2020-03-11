#' svm_test allows assessing the final DEGs through a machine learning step by using SVM with a test dataset.
#'
#' svm_test allows assessing the final DEGs through a machine learning step by using SVM with a test dataset. An optimization of C and G hiperparameters is done at the start of the process.
#'
#' @param test The test parameter is an expression matrix or data.frame that contains the test dataset with the genes in the columns and the samples in the rows.
#' @param labelsTest A vector or factor that contains the test labels for each of the samples in the test object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param bestParameters Best values for C and gamma parameters selected during the training phase.
#' @return A list that contains four objects. The confusion matrix, the accuracy, the sensitibity and the specificity for each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' testMatrix <- t(DEGsMatrix)[c(5,10),]
#' testLabels <- labels[c(5,10)]
#' #previously trained svm model
#' bestParameters <- c(cost = Rsvm_sb$bestTune$C, gamma = Rsvm_sb$bestTune$sigma)
#' svm_test(testMatrix,testLabels,rownames(DEGsMatrix)[1:10], bestParameters)

svm_test <-function(test,labelsTest,vars_selected,bestParameters){

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

  for(i in seq_len(dim(train)[2])){

    cat(paste("Testing with ", i," variables...\n",sep=""))
    svm_model<-svm(train[,seq(i)],labelsTrain,kernel='radial',
                   cost=bestParameters$C,gamma=bestParameters$gamma,probability=TRUE)
    predicts<-predict(svm_model,test[,seq(i)],probability=TRUE)

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
