#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest.
#'
#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest. Furthermore, the ranking is returned and can be used as input of the parameter vars_selected in the machine learning functions.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each samples in data parameter.
#' @param vars_selected The genes selected to use in the feature selection process. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes.
#' @param mode The algorithm used to calculate the genes ranking. The possibilities are two: mrmr, rf and da.
#' @param disease The name of a disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @param subdiseases Vector with the name of the particular subdiseases from disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @return A vector that contains the ranking of genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' featureRanking <- featureSelection(t(DEGsMatrix),labels,rownames(DEGsMatrix),mode='mrmr')
#' featureRanking <- featureSelection(t(DEGsMatrix),labels,rownames(DEGsMatrix),mode='daRed',disease='cancer')
#' featureRanking <- featureSelection(t(DEGsMatrix),labels,rownames(DEGsMatrix),mode='daRed',disease='cancer',subdiseases = c('colorectal cancer','breast cancer'))


featureSelection <-function(data,labels,vars_selected,mode="mrmr",disease="",subdiseases=c()){

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

  }else if(mode == "da" || mode == 'daRed'){
    
    if(disease == ""){
      stop("Please, indicate a disease name to acquire the Disease Association Score and Feature selection.")
    }
    
    cat("Calculating ranking of biological relevant genes by using DA implementation...\n")
    
    relatedDiseases <- DEGsToDiseases(vars_selected, size = 100)
    
    # Check if it is not a multiclass problem
    if (length(subdiseases) == 0){
      overallRanking <- c()
      for(i in seq(length(relatedDiseases))){
        
        if(is.na(grep(disease,relatedDiseases[[i]]$summary[,1])[1]))
          overallScore <- 0.0
        else
          overallScore <- relatedDiseases[[i]]$summary[grep(disease,relatedDiseases[[i]]$summary[,1])[1],2]

        overallRanking <- c(overallRanking,overallScore)
      }
    }
    else{
      overallRankingDisease <- list()
      # Get DA score for each sub disease
      for( subdisease in subdiseases){
        overallRankingDisease[[subdisease]] = c()
        for(i in seq(length(relatedDiseases))){
          if(is.na(grep(subdisease,relatedDiseases[[i]]$summary[,1])[1])){
            overallScore <- 0.0
          }else{
            overallScore <- relatedDiseases[[i]]$summary[grep(subdisease,relatedDiseases[[i]]$summary[,1])[1],2]
          }
          overallRankingDisease[[subdisease]] <- c(overallRankingDisease[[subdisease]],overallScore)
        }
        names(overallRankingDisease[[subdisease]]) <- names(relatedDiseases)
      }
      
      # Final DA score is the difference in absolute value
      overallRanking <- rep(0,length(overallRankingDisease[[1]]))
      for (i in seq(length(overallRankingDisease))){
        overallRanking <- abs(overallRanking-as.numeric(overallRankingDisease[[i]]))
      }
    }
    
    names(overallRanking) <- names(relatedDiseases)
    overallRanking <- sort(overallRanking,decreasing = TRUE)
    
    cat("Disease Association ranking: ")
    cat(names(overallRanking))
    cat("\n")
    
    if (mode == 'da') return(overallRanking)
    
    if (length(subdiseases) == 0){
      evidences <- DEGsEvidences(names(overallRanking), disease, size=100)
      cat("Calculating redundances between found evidences...\n")
      redundances <- evidencesToRedundance(evidences)
    }else{
      evidences <- list()
      act.redundances <- list()
      redundances <- matrix(0,ncol=length(overallRanking),nrow=length(overallRanking))

      for (subdisease in subdiseases){
        evidences[[subdisease]] <- DEGsEvidences(names(overallRanking),disease,subdisease,size=10)
        cat(paste("Calculating redundances between found evidences for subdisease",disease,"...\n"))
        act.redundances[[subdisease]] <- evidencesToRedundance(evidences[[subdisease]])
        redundances <- abs(redundances - act.redundances[[subdisease]])
      }
      colnames(redundances) = colnames(act.redundances[[subdisease[1]]])
      rownames(redundances) = rownames(act.redundances[[subdisease[1]]])
    }
    
    
    cat("Calculating genes scores...\n")
    genes <- names(overallRanking)
    
    # Output: list of selected genes
    selected.genes <- list()
    # Select a gene with maximun score
    selected.genes[[genes[which(overallRanking == max(overallRanking))[1]]]] = max(overallRanking)
    
    # Iter over all genes
    for (i in seq(length(overallRanking)-1)){
      # actual max score
      max <- -1000
      
      # act.genes contains genes to select (genes that are not already selected)
      act.genes <- genes[ ! genes %in% names(selected.genes) ]
      
      # Iter in genes to select
      for ( gen in act.genes){
        # Gene relevance is it's score
        rel <- as.numeric(overallRanking[[gen]])
        # Redundance begin as 0
        red <- 0
        
        # Calculate redundance between actual gene and selected genes
        for ( sel in names(selected.genes)){
          red <- red + redundances[gen,sel]
        }
        # Score = relevance - relevance * redundance/ num of selected genes
        score <- rel - red/length(selected.genes) * rel
        
        # If actual score is max. keep actual gene
        if (score > max){
          max <- score
          act.selected <- gen
        }
      }
      # Save gene with maximun score
      selected.genes[[act.selected]] <- max
    }
    return(selected.genes)
    
  }else{
    stop("The mode is unrecognized, please use mrmr, rf, da or daRed.")
  }
}


