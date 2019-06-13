#' rf_CV allows assessing the final DEGs through a machine learning step by using Random Forest in a cross validation process.
#'
#' rf_CV allows assessing the final DEGs through a machine learning step by using Random Forest in a cross validation process. This function applies a cross validation of n folds with representation of all classes in each fold. The 80\% of the data are used for training and the 20\% for test.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param numFold The number of folds to carry out in the cross validation process.
#' @return A list that contains four objects. The confusion matrix for each fold, the accuracy, the sensitibity and the specificity for each fold and each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' rf_CV(t(DEGsMatrix),labels,rownames(DEGsMatrix),2)

rf_CV<-function(data,labels,vars_selected,numFold=10){


  if(!is.data.frame(data) && !is.matrix(data)){

    stop("The data argument must be a dataframe or a matrix.")

  }
  if(dim(data)[1] != length(labels)){

    stop("The length of the rows of the argument data must be the same than the length of the lables. Please, ensures that the rows are the samples and the columns are the variables.")

  }

  if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
  if(is.character(labels)){ labels <- as.factor(labels) }

  if(numFold%%1!=0 || numFold == 0){

    stop("The numFold argument must be integer and greater than 0.")

  }

  data <- as.data.frame(apply(data,2,as.double))
  data <- data[,vars_selected]

  data = vapply(data, function(x){ 
    max = max(x)
    min = min(x)
    x = ((x-min)/(max-min))*2-1}, double(nrow(data)))
  
  data <- as.data.frame(data)
  
  acc_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  sens_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  spec_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])

  cfMatList  <- list()

  for(i in seq_len(numFold)){

    cat(paste("Training fold ", i,"...\n",sep=""))
    trainingDataset <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 0)),
                                colnames(data))
    testDataset <- setNames(data.frame(matrix(ncol = ncol(data), nrow = 0)),
                            colnames(data))
    labelsTrain <- factor(0L)
    labelsTest <- factor(0L)

    for(class in names(table(labels))){

      classPos <- which(labels == class)
      classPos <- sample(classPos)
      trainingPos <- round(length(classPos)*0.8)
      testPos <- round(length(classPos)*0.2)
      trainingDataset <- rbind(trainingDataset,data[classPos[seq_len(trainingPos)],])
      testDataset <- rbind(testDataset,data[classPos[(trainingPos+1):(trainingPos+testPos)],])
      labelsTrain <- unlist(list(labelsTrain, labels[classPos[seq_len(trainingPos)]]))
      labelsTest <- unlist(list(labelsTest, labels[classPos[(trainingPos+1):(trainingPos+testPos)]]))

    }

    labelsTrain <- factor(labelsTrain[-1])
    labelsTest <- factor(labelsTest[-1])

    for(j in 2:length(vars_selected)){

      rf_mod = randomForest(x = trainingDataset[,seq(j)], y = labelsTrain, ntree = 100)
      predicts <- predict(rf_mod , testDataset[,seq(j)])
      cfMatList[[i]] <- confusionMatrix(predicts,labelsTest)
      acc_cv[i,j]<-confusionMatrix(predicts,labelsTest)$overall[[1]]
      sens_cv[i,j]<-confusionMatrix(predicts,labelsTest)$byClass[[1]]
      spec_cv[i,j]<-confusionMatrix(predicts,labelsTest)$byClass[[2]]

    }
  }

  rownames(acc_cv) <- paste("Fold",seq(numFold),sep = "")
  colnames(acc_cv) <- vars_selected
  rownames(sens_cv) <- paste("Fold",seq(numFold),sep = "")
  colnames(sens_cv) <- vars_selected
  rownames(spec_cv) <- paste("Fold",seq(numFold),sep = "")
  colnames(spec_cv) <- vars_selected

  cat("Classification done successfully!\n")
  results_cv <- list(cfMatList,acc_cv,sens_cv,spec_cv)
  names(results_cv) <- c("cfMats","accMatrix","sensMatrix","specMatrix")
  invisible(results_cv)

}
