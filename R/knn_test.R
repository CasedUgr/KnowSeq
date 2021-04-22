#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset.
#'
#' knn_test allows assessing the final DEGs through a machine learning step by using k-NN with a test dataset. An optimization of the k neighbours is done at the start of the process.
#'
#' @param train The train parameter is an expression matrix or data.frame that contains the train dataset with the genes in the columns and the samples in the rows.
#' @param labelsTrain A vector or factor that contains the train labels for each of the samples in the train object.
#' @param test The test parameter is an expression matrix or data.frame that contains the test dataset with the genes in the columns and the samples in the rows.
#' @param labelsTest A vector or factor that contains the test labels for each of the samples in the test object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{DEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param bestK Best K selected during the training phase.
#' @return A list that contains four objects. The confusion matrix, the accuracy, the sensitibity and the specificity for each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' trainingMatrix <- t(DEGsMatrix)[c(1:4,6:9),]
#' trainingLabels <- labels[c(1:4,6:9)]
#' testMatrix <- t(DEGsMatrix)[c(5,10),]
#' testLabels <- labels[c(5,10)]
#' bestK <- 3 # the one that has been selected
#' results_test_knn <- knn_test(trainingMatrix, trainingLabels, testMatrix, testLabels, rownames(DEGsMatrix)[1:10], bestK)

knn_test <-function(train,labelsTrain,test,labelsTest,vars_selected, bestK){

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
    max <- max(x)
    min <- min(x)
    if(max >  min){
      x <- ((x - min) / (max - min)) * 2 - 1
    }
    else{
      x
    }}, double(nrow(train)))
  
  train <- as.data.frame(train)
  
  test = vapply(test, function(x){ 
    max <- max(x)
    min <- min(x)
    if(max >  min){
      x <- ((x - min) / (max - min)) * 2 - 1
    }
    else{
      x
    }}, double(nrow(test)))
  
  test <- as.data.frame(test)

  accVector <- double()
  sensVector <- double()
  specVector <- double()
  f1Vector <- double()
  cfMatList  <- list()
  predictsVector <- list()
  
  # Firstly with one variable
  cat(paste("Testing with ", 1," variables...\n",sep=""))
  knn_mod = knn3(x = train[, 1, drop=FALSE], y = labelsTrain, k = bestK)
  predictScores <- predict(knn_mod, test[, 1, drop=FALSE], type = "prob")
  predicts <- predict(knn_mod, test[, 1, drop=FALSE], type = "class")
  
  cfMat<-confusionMatrix(predicts,labelsTest)
  if (length(levels(labelsTrain))==2){
    sens <- cfMat$byClass[[1]]
    spec <- cfMat$byClass[[2]]
    f1 <- cfMat$byClass[[7]]
  } else{
    sens <- mean(cfMat$byClass[,1])
    spec <- mean(cfMat$byClass[,2])
    f1 <- mean(cfMat$byClass[,7])
  }

  cfMatList[[1]] <- cfMat
  accVector[1] <- cfMat$overall[[1]]
  sensVector[1] <- sens
  specVector[1] <- spec
  f1Vector[1] <- f1
  predictsVector[[1]] <- predictScores
  if(is.na(f1Vector[1])) f1Vector[i] <- 0
  
  for(i in c(2:dim(test)[2])){

    cat(paste("Testing with ", i," variables...\n",sep=""))
    knn_mod = knn3(x = train[,seq(i)], y = labelsTrain, k = bestK)
    predictScores <- predict(knn_mod, test[,seq(i)], type = "prob")
    predicts <- predict(knn_mod, test[,seq(i)], type = "class")
    
    cfMat<-confusionMatrix(predicts,labelsTest)
    
    if (length(levels(labelsTrain))==2){
      sens <- cfMat$byClass[[1]]
      spec <- cfMat$byClass[[2]]
      f1 <- cfMat$byClass[[7]]
    } else{
      sens <- mean(cfMat$byClass[,1])
      spec <- mean(cfMat$byClass[,2])
      f1 <- mean(cfMat$byClass[,7])
    }
    
    cfMatList[[i]] <- cfMat
    accVector[i] <- cfMat$overall[[1]]
    sensVector[i] <- sens
    specVector[i] <- spec
    f1Vector[i] <- f1
    predictsVector[[i]] <- predictScores
    if(is.na(f1Vector[i])) f1Vector[i] <- 0
  }

  cat("Classification done successfully!\n")
  names(accVector) <- vars_selected
  names(sensVector) <- vars_selected
  names(specVector) <- vars_selected
  names(f1Vector) <- vars_selected

  results <- list(cfMatList,accVector,sensVector,specVector,f1Vector,predictsVector)
  names(results) <- c("cfMats","accVector","sensVector","specVector","f1Vector","predictions")
  invisible(results)

}
