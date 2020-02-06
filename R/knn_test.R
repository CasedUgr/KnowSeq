#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset.
#'
#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset. An optimization of the k neighbours is done at the start of the process.
#'
#' @param train The train parameter is an expression matrix or data.frame that contains the training dataset with the genes in the columns and the samples in the rows.
#' @param labelsTrain A vector or factor that contains the training labels for each of the samples in the train object.
#' @param test The test parameter is an expression matrix or data.frame that contains the test dataset with the genes in the columns and the samples in the rows.
#' @param labelsTest A vector or factor that contains the test labels for each of the samples in the test object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @return A list that contains four objects. The confusion matrix, the accuracy, the sensitibity and the specificity for each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' trainingMatrix <- t(DEGsMatrix)[c(1:4,6:9),]
#' trainingLabels <- labels[c(1:4,6:9)]
#' testMatrix <- t(DEGsMatrix)[c(5,10),]
#' testLabels <- labels[c(5,10)]
#'
#' results_test_knn <- knn_test(trainingMatrix, trainingLabels, testMatrix, testLabels, rownames(DEGsMatrix)[1:10])

knn_test <-function(train,labelsTrain,test,labelsTest,vars_selected){

  if(!is.data.frame(train) && !is.matrix(train)){

    stop("The train argument must be a dataframe or a matrix.")

  }

  if(dim(train)[1] != length(labelsTrain)){

    stop("The length of the rows of the argument train must be the same than the length of the lablesTrain. Please, ensures that the rows are the samples and the columns are the variables.")

  }

  if(!is.character(labelsTrain)  && !is.factor(labelsTrain)){stop("The class of the labelsTrain parameter must be character vector or factor.")}
  if(is.character(labelsTrain)){ labelsTrain <- as.factor(labelsTrain) }

  if(!is.character(labelsTest)  && !is.factor(labelsTest)){stop("The class of the labelsTest parameter must be character vector or factor.")}
  if(is.character(labelsTest)){ labelsTest <- as.factor(labelsTest) }

  if(!is.data.frame(test) && !is.matrix(test)){

    stop("The test argument must be a dataframe or a matrix.")

  }

  if(dim(test)[1] != length(labelsTest)){

    stop("The length of the rows of the argument test must be the same than the length of the lablesTest. Please, ensures that the rows are the samples and the columns are the variables.")

  }

  train <- as.data.frame(apply(train,2,as.double))
  train <- train[,vars_selected]
  test <- as.data.frame(apply(test,2,as.double))
  test <- test[,vars_selected]

  train = vapply(train, function(x){ 
    max = max(x)
    min = min(x)
    x = ((x-min)/(max-min))*2-1}, double(nrow(train)))
  
  train <- as.data.frame(train)
  
  test = vapply(test, function(x){ 
    max = max(x)
    min = min(x)
    x = ((x-min)/(max-min))*2-1}, double(nrow(test)))
  
  test <- as.data.frame(test)

  fitControl <- trainControl(method = "repeatedcv", number = 10,repeats=3)
  cat("Tuning the optimal K...\n")
  K_sb <- train(train, labelsTrain, method = "knn",trControl = fitControl,preProcess = c("center", "scale"),tuneLength = 10)

  bestK = K_sb$bestTune

  accVector <- double()
  sensVector <- double()
  specVector <- double()
  cfMatList  <- list()

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
