#' svm_CV allows assessing the final DEGs through a machine learning step by using svm in a cross validation process.
#'
#' svm_CV allows assessing the final DEGs through a machine learning step by using svm in a cross validation process. This function applies a cross validation of n folds with representation of all classes in each fold. The 80\% of the data are used for training and the 20\% for test. An optimization of C and G hiperparameters is done at the start of the process.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each of the samples in the data object.
#' @param vars_selected The genes selected to classify by using them. It can be the final DEGs extracted with the function \code{\link{DEGsExtraction}} or a custom vector of genes. Furthermore, the ranking achieved by \code{\link{featureSelection}} function can be used as input of this parameter.
#' @param numFold The number of folds to carry out in the cross validation process.
#' @return A list that contains five objects. The confusion matrix for each fold, the accuracy, the sensitibity and the specificity for each fold and each genes, and a vector with the best parameters found for the SVM algorithm after tuning.
#' @examples
#' dir <- system.file("extdata", package = "KnowSeq")
#' load(paste(dir, "/expressionExample.RData", sep = ""))
#'
#' svm_CV(t(DEGsMatrix), labels, rownames(DEGsMatrix), 2)

svm_CV <- function(data, labels, vars_selected, numFold = 10) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("The data argument must be a dataframe or a matrix.")
  }
  if (dim(data)[1] != length(labels)) {
    stop("The length of the rows of the argument data must be the same than the length of the lables. Please, ensures that the rows are the samples and the columns are the variables.")
  }

  if (!is.character(labels) && !is.factor(labels)) {
    stop("The class of the labels parameter must be character vector or factor.")
  }
  if (is.character(labels)) {
    labels <- as.factor(labels)
  }

  if (numFold %% 1 != 0 || numFold == 0) {
    stop("The numFold argument must be integer and greater than 0.")
  }

  data <- as.data.frame(apply(data, 2, as.double))
  data <- data[, vars_selected]

  data <- vapply(data, function(x) {
    max <- max(x)
    min <- min(x)
    x <- ((x - min) / (max - min)) * 2 - 1
  }, double(nrow(data)))

  data <- as.data.frame(data)

  fitControl <- caret::trainControl(method = "cv", number = 10)
  cat("Tuning the optimal C and G...\n")

  grid_radial <- expand.grid(
    sigma = c(
      0, 0.01, 0.02, 0.025, 0.03, 0.04,
      0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.25, 0.5, 0.75, 0.9
    ),
    C = c(
      0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
      1, 1.5, 2, 5
    )
  )
  
  dataForTunning <- cbind(data, labels)

  Rsvm_sb <- train(labels ~ ., data = dataForTunning, type = "C-svc", method = "svmRadial", preProc = c("center", "scale"), trControl = fitControl, tuneGrid = grid_radial)
  
  bestParameters <- c(cost = Rsvm_sb$bestTune$C, gamma = Rsvm_sb$bestTune$sigma)
  cat(paste("Optimal cost:", bestParameters[1], "\n"))
  cat(paste("Optimal gamma:", bestParameters[2], "\n"))
  
  acc_cv <- matrix(0L, nrow = numFold, ncol = dim(data)[2])
  sens_cv <- matrix(0L, nrow = numFold, ncol = dim(data)[2])
  spec_cv <- matrix(0L, nrow = numFold, ncol = dim(data)[2])
  cfMatList <- list()
  # compute size of val fold
  lengthValFold <- dim(data)[1]/numFold
  
  # reorder the data matrix in order to have more
  # balanced folds
  positions <- rep(seq_len(dim(data)[1]))
  randomPositions <- sample(positions)
  data <- data[randomPositions,]
  labels <- labels[randomPositions]
  
  for (i in seq_len(numFold)) {
    cat(paste("Training fold ", i, "...\n", sep = ""))
    
    # obtain validation and training folds
    valFold <- seq(round((i-1)*lengthValFold + 1 ), round(i*lengthValFold))
    trainDataCV <- setdiff(seq_len(dim(data)[1]), valFold)
    testDataset<- data[valFold,]
    trainingDataset <- data[trainDataCV,]
    labelsTrain <- labels[trainDataCV]
    labelsTest <- labels[valFold]

    for (j in seq_len(length(vars_selected))) {
      svm_model <- svm(trainingDataset[, seq(j)], labelsTrain,
        kernel = "radial",
        cost = Rsvm_sb$bestTune$C, gamma = Rsvm_sb$bestTune$sigma, probability = TRUE
      )
      predicts <- predict(svm_model, testDataset[, seq(j)], probability = TRUE)

      cfMatList[[i]] <- confusionMatrix(predicts, labelsTest)
      acc_cv[i, j] <- confusionMatrix(predicts, labelsTest)$overall[[1]]
      sens_cv[i, j] <- confusionMatrix(predicts, labelsTest)$byClass[[1]]
      spec_cv[i, j] <- confusionMatrix(predicts, labelsTest)$byClass[[2]]
      
      if(is.na(sens_cv[i,j])) sens_cv[i,j] <- 0
      if(is.na(spec_cv[i,j])) spec_cv[i,j] <- 0
    }
  }
  rownames(acc_cv) <- paste("Fold", seq(numFold), sep = "")
  colnames(acc_cv) <- vars_selected
  rownames(sens_cv) <- paste("Fold", seq(numFold), sep = "")
  colnames(sens_cv) <- vars_selected
  rownames(spec_cv) <- paste("Fold", seq(numFold), sep = "")
  colnames(spec_cv) <- vars_selected

  cat("Classification done successfully!\n")
  results_cv <- list(cfMatList, acc_cv, sens_cv, spec_cv, bestParameters)
  names(results_cv) <- c("cfMats", "accMatrix", "sensMatrix", "specMatrix", "bestParameters")
  invisible(results_cv)
}
