#' rf_trn allows assessing the final DEGs through a machine learning step by using Random Forest in a cross validation process.
#'
#' rf_trn allows assessing the final DEGs through a machine learning step by using Random Forest in a cross validation process. This function applies a cross validation of n folds with representation of all classes in each fold. The 80\% of the data are used for training and the 20\% for test.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{DEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param numFold The number of folds to carry out in the cross validation process.
#' @return A list that contains four objects. The confusion matrix for each fold, the accuracy, the sensitibity and the specificity for each fold and each genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' rf_trn(t(DEGsMatrix)[,1:10],labels,rownames(DEGsMatrix)[1:10],2)

rf_trn <- function(data,labels,vars_selected,numFold=10){


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
  
  fitControl <- trainControl(method = "cv", number = 10)
  cat("Tuning the optimal mtry...\n")
  
  dataForTunning <- cbind(data, labels)
  tunegrid <- expand.grid(.mtry = (1:30)) 
 
  rf_gridsearch <- train(labels ~ ., 
                         data = dataForTunning,
                         method = 'rf',
                         metric = 'Accuracy',
                         preProc = c("center", "scale"),
                         ntree=1000,
                         tuneGrid = tunegrid)
  mtryTune <- rf_gridsearch$bestTune$mtry
  
  acc_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  sens_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  spec_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  f1_cv<-matrix(0L,nrow = numFold,ncol = dim(data)[2])
  
  cfMatList  <- list()
  # compute size of val fold
  lengthValFold <- dim(data)[1]/numFold
  
  # reorder the data matrix in order to have more
  # balanced folds
  positions <- rep(seq_len(dim(data)[1]))
  randomPositions <- sample(positions)
  data <- data[randomPositions,]
  labels <- labels[randomPositions]
  
  for(i in seq_len(numFold)){

    cat(paste("Training fold ", i,"...\n",sep=""))
    
    # obtain validation and training folds
    valFold <- seq(round((i-1)*lengthValFold + 1 ), round(i*lengthValFold))
    trainDataCV <- setdiff(seq_len(dim(data)[1]), valFold)
    testDataset<- data[valFold,]
    trainingDataset <- data[trainDataCV,]
    labelsTrain <- labels[trainDataCV]
    labelsTest <- labels[valFold]
    colNames <- colnames(trainingDataset)
    
    # first iteration is performed outside of the foor lopp
    # in order to avoid having a if inside
    columns <- c(colNames[1])
    tr_ctr <- trainControl(method="none")
    dataForTrt <- data.frame(cbind(subset(trainingDataset, select=columns),labelsTrain))
    colnames(dataForTrt)[seq(1)] <- columns
    rf_mod <- train(labelsTrain ~ ., 
          data = dataForTrt,
          method = 'rf',
          metric = 'Accuracy',
          preProc = c("center", "scale"),
          ntree=1000,
          tuneGrid = data.frame(.mtry= mtryTune))
    
    unkX <- subset(testDataset, select=columns)
    predicts <- extractPrediction(list(my_rf=rf_mod), testX = subset(testDataset, select=columns), unkX = unkX,
                                  unkOnly = !is.null(unkX) & !is.null(subset(testDataset, select=columns)))
    
    predicts <- predicts$pred

    cfMatList[[i]] <- confusionMatrix(predicts,labelsTest)
    acc_cv[i,1]<-cfMatList[[i]]$overall[[1]]
    
    if (length(levels(labelsTrain))==2){
      sens <- cfMatList[[i]]$byClass[[1]]
      spec <- cfMatList[[i]]$byClass[[2]]
      f1 <- cfMatList[[i]]$byClass[[7]]
    } else{
      sens <- mean(cfMatList[[i]]$byClass[,1])
      spec <- mean(cfMatList[[i]]$byClass[,2])
      f1 <- mean(cfMatList[[i]]$byClass[,7])
    }
    
    sens_cv[i,1]<-sens
    spec_cv[i,1]<-spec
    f1_cv[i,1]<-f1
    
    if(is.na(sens_cv[i,1])) sens_cv[i,1] <- 0
    if(is.na(spec_cv[i,1])) spec_cv[i,1] <- 0
    if(is.na(f1_cv[i,1])) f1_cv[i,1] <- 0
    
    for(j in 2:length(vars_selected)){

      columns <- c(colNames[seq(j)])
      tr_ctr <- trainControl(method="none")
      dataForTrt <- data.frame(cbind(subset(trainingDataset, select=columns),labelsTrain))
      colnames(dataForTrt)[seq(j)] <- columns
      rf_mod <- train(labelsTrain ~ ., 
                      data = dataForTrt,
                      method = 'rf',
                      metric = 'Accuracy',
                      preProc = c("center", "scale"),
                      ntree=1000,
                      tuneGrid = data.frame(.mtry= mtryTune))
      
      unkX <- subset(testDataset, select=columns)
      predicts <- extractPrediction(list(my_rf=rf_mod), testX = subset(testDataset, select=columns), unkX = unkX,
                                    unkOnly = !is.null(unkX) & !is.null(subset(testDataset, select=columns)))
      
      predicts <- predicts$pred
      cfMatList[[i]] <- confusionMatrix(predicts,labelsTest)
      acc_cv[i,j]<-cfMatList[[i]]$overall[[1]]
      
      if (length(levels(labelsTrain))==2){
        sens <- cfMatList[[i]]$byClass[[1]]
        spec <- cfMatList[[i]]$byClass[[2]]
        f1 <- cfMatList[[i]]$byClass[[7]]
      } else{
        sens <- mean(cfMatList[[i]]$byClass[,1])
        spec <- mean(cfMatList[[i]]$byClass[,2])
        f1 <- mean(cfMatList[[i]]$byClass[,7])
      }
      
      sens_cv[i,j]<-sens
      spec_cv[i,j]<-spec
      f1_cv[i,j]<-f1
      
      if(is.na(sens_cv[i,j])) sens_cv[i,j] <- 0
      if(is.na(spec_cv[i,j])) spec_cv[i,j] <- 0
      if(is.na(f1_cv[i,j])) f1_cv[i,j] <- 0

    }
  }

  meanAcc <- colMeans(acc_cv)
  names(meanAcc) <- colnames(acc_cv)
  sdAcc <- apply(acc_cv, 2, sd)
  accuracyInfo <- list(meanAcc, sdAcc)
  names(accuracyInfo) <- c("meanAccuracy","standardDeviation")
  
  
  meanSens <- colMeans(sens_cv)
  names(meanSens) <- colnames(sens_cv)
  sdSens <- apply(sens_cv, 2, sd)
  sensitivityInfo <- list(meanSens, sdSens)
  names(sensitivityInfo) <- c("meanSensitivity","standardDeviation")
  
  
  meanSpec <- colMeans(spec_cv)
  names(meanSpec) <- colnames(spec_cv)
  sdSpec <- apply(spec_cv, 2, sd)
  specificityInfo <- list(meanSpec, sdSpec)
  names(specificityInfo) <- c("meanSpecificity","standardDeviation")
  
  
  meanF1 <- colMeans(f1_cv)
  names(meanF1) <- colnames(f1_cv)
  sdF1 <- apply(f1_cv, 2, sd)
  F1Info <- list(meanF1, sdF1)
  names(F1Info) <- c("meanF1","standardDeviation")
  
  cat("Classification done successfully!\n")
  results_cv <- list(cfMatList,accuracyInfo,sensitivityInfo,specificityInfo,F1Info,mtryTune)
  names(results_cv) <- c("cfMats","accuracyInfo","sensitivityInfo","specificityInfo","F1Info","bestParameters")
  invisible(results_cv)

}
