#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest.
#'
#' featureSelection function calculates the optimal order of DEGs to achieve the best result in the posterior machine learning process by using mRMR algorithm or Random Forest. Furthermore, the ranking is returned and can be used as input of the parameter vars_selected in the machine learning functions.
#'
#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each samples in data parameter.
#' @param vars_selected The genes selected to use in the feature selection process. It can be the final DEGs extracted with the function \code{\link{DEGsExtraction}} or a custom vector of genes.
#' @param mode The algorithm used to calculate the genes ranking. The possibilities are two: mrmr, rf and da.
#' @param disease The name of a disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @param subdiseases Vector with the name of the particular subdiseases from disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @param maxGenes Integer that indicated the maximun number of genes to be returned. 
#' @return A vector that contains the ranking of genes.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#' featureRanking <- featureSelection(t(DEGsMatrix),labels,rownames(DEGsMatrix),mode='mrmr')

featureSelection <-function(data,labels,vars_selected,mode="mrmr",disease="",subdiseases=c(),maxGenes=ncol(data)){
  
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
    if(mode == "da"){
      cat("Calculating ranking of biological relevant genes by using DA implementation...\n")
    }else if(mode == 'daRed'){
      cat("Calculating ranking of biological relevant genes by using DA-Red implementation...\n")
    }
    
    
    if (length(subdiseases)==0){
      overallRanking <- rep(0,length(vars_selected))
      names(overallRanking) <- vars_selected
      
      disease_ <-  str_replace_all(disease,' ','-')
      r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease_,"&size=1&filter=disease",sep = ""))
      respon <- content(r_Ensembl)
      
      if ( 'size' %in% names(respon) && respon$size == 0) stop('Disease not found')
      
      disease.id <- respon$data[[1]]$id
      url  <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=",disease.id,"&size=10000",sep='')
      response <- GET(url)
      response <- content(response)
      found.symbols <- unlist(list.map(response$data,target$gene_info$symbol))
      scores <- unlist(list.map(response$data,association_score$overall))
      
      keep <- which(!duplicated(found.symbols))
      found.symbols[keep]
      scores[keep]

      found.index <- which(found.symbols %in% vars_selected)
      overallRanking[found.symbols[found.index]] <- scores[found.index]
    }
    else{
      overallRanking <- rep(0,length(vars_selected))
      names(overallRanking) <- vars_selected
      
      for (i in seq(length(subdiseases))){
        subdisease <- subdiseases[i]
        subdisease_ <- str_replace_all(subdisease)
        r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",subdisease_,"&size=1&filter=disease",sep = ""))
        respon <- content(r_Ensembl)

        if ( 'size' %in% names(respon) && respon$size == 0) stop(paste('Subdisease',subdisease,'not found'))

        disease.id <- respon$data[[1]]$id
        url  <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=",disease.id,"&size=10000",sep='')
        response <- GET(url)
        response <- content(response)
        found.symbols <- unlist(list.map(response$data,target$gene_info$symbol))
        scores <- unlist(list.map(response$data,association_score$overall))

        keep <- which(!duplicated(found.symbols))
        found.symbols <- found.symbols[keep]
        scores <- scores[keep]
        
        found.index <- which(found.symbols %in% vars_selected)
        if (i==1) overallRanking[found.symbols[found.index]] <- as.vector(scores[found.index])
        else overallRanking[found.symbols[found.index]] <- rowMeans(cbind(as.vector(overallRanking[found.symbols[found.index]]),as.vector(scores[found.index])))
      }
    }

    overallRanking <- sort(overallRanking,decreasing = TRUE)
    cat("Disease Association ranking: ")
    cat(names(overallRanking))
    cat("\n")
    
    if (mode == 'da') return(overallRanking)
    
    if (length(subdiseases) == 0){
      
      evidences <- DEGsEvidences(names(overallRanking), disease, size=100)
      
    }else{
      evidences <- list()
      
      for (subdisease in subdiseases){
        
        evidences[[subdisease]] <- DEGsEvidences(names(overallRanking),subdisease,size=100)
        
      }
      remove <- list('global'=c())
      for ( gen in names(names(overallRanking))){
        all = TRUE
        for (subdisease in subdiseases){
          if (! subdisease %in% names(evidences)) evidences[[subdisease]] <- c()
          if (is.character(evidences[[subdisease]][[gen]])) remove[[subdisease]] <- c(remove[[subdisease]],gen)
          else all = FALSE
        }
        if ( all ) remove[['global']] <- c(remove[['global']],gen)
      }
      
      # Join subdiseases evidences 
      aux <- list()
      for (subdisease in subdiseases){
        evidences[[subdisease]] <- evidences[[subdisease]][! names(evidences[[subdisease]]) %in% remove[[subdisease]]]
        aux <- list(aux,evidences[[subdisease]])
      }
      
      final.evidences <- evidences[[subdiseases[1]]]
      for (gen in names(final.evidences)){
        if ( is.list(evidences[[subdiseases[2]]][[gen]])){
          for (type in names(evidences[[subdiseases[2]]][[gen]])){
            if (type %in% names(final.evidences[[gen]])){
              final.evidences[[gen]][[type]] <- append(final.evidences[[gen]][[type]],evidences[[subdiseases[2]]][[gen]][[type]])
            }
            else{
              if (is.list(final.evidences[[gen]])) final.evidences[[gen]]  <-  list()
              final.evidences[[gen]][[type]] <- evidences[[subdiseases[2]]][[gen]][[type]]
            }
          }
        }
      }
      evidences <- final.evidences
      for (gen in remove[['global']]){
        evidences[[gen]] <- 'Not evidences found'
      }
    }
    
    cat("Calculating genes scores...\n")
    genes <- names(overallRanking)
    
    # Output: list of selected genes
    selected.genes <- list()
    # Select a gene with maximun score
    selected.genes[[names(overallRanking)[1]]] = overallRanking[1]
    
    # Create empty redundance matrix
    redundances <- matrix(-1,ncol=maxGenes,nrow=length(overallRanking))
    rownames(redundances) = names(overallRanking)
    colnames(redundances) <- rep('',maxGenes)
    colnames(redundances)[1] <- names(selected.genes[1])
    
    # Iter over all genes
    for (i in seq(maxGenes-1)){
      # actual max score
      max.score <- -1000
      
      # act.genes contains genes to select (genes that are not already selected)
      
      act.genes <- genes[ ! genes %in% names(selected.genes) ]
      
      # Iter in genes to select
      for ( gen in act.genes){
        # Gene relevance is it's score
        rel <- as.numeric(overallRanking[[gen]])
        # Redundance begin as 0
        red <- 0
        
        
        if (rel < max.score) break
        
        # Calculate redundance between actual gene and selected genes
        for ( sel in names(selected.genes)){
          # Calculate and save redudance 
          if (redundances[gen,sel] == -1){
            redundances[gen,sel] = 0
            if ( is.list(evidences[[gen]]) && is.list(evidences[[sel]])){
              gen.nevs <- 0
              for (type in names(evidences[[gen]])){
                if (type %in% names(evidences[[gen]])){
                  # type.total contains found coincidences for actual type of evidences
                  type.total  <- 0
                  ncol <- length(evidences[[gen]][[type]][[1]]$evidence)
                  # Iter on gen evidences
                  for (row1 in evidences[[gen]][[type]]){
                    gen.nevs <- gen.nevs + 1
                    # Boolean matrix that contains coincidences between gen and sel
                    
                    act.total.matrix <- matrix(0,nrow=length(evidences[[sel]][[type]]),ncol=ncol)
                    # Iter on gen evidences
                    for (row2 in evidences[[sel]][[type]]){
                      # Add row to act.total.matrix with boolean values
                      act.total.matrix <- rbind(act.total.matrix, unlist(row1$evidence) == unlist(row2$evidence))
                      # This row fully coincide with actual evidencie, so we stop searching
                      if (rowSums(tail(act.total.matrix,1) == ncol))  break
                    }
                    # If there are any row that fully coincide add 1 (this evidences is fully contained in sel data)
                    if ( any(rowSums(act.total.matrix) == ncol)) type.total = type.total + 1
                    else{
                      # Check which row is the most similar to the actual evidence
                      # Add the percentage of coincidence ( num. coincidences / ncol)
                      act.total.field <- colSums(act.total.matrix)
                      type.total = type.total + length(which(colSums(act.total.matrix) >= 1))/ncol
                    }
                  }
                  # Add score of type evidences
                  redundances[gen,sel] = redundances[gen,sel] + type.total
                }
              }
              # Normalize score dividing by the number of evidences
              redundances[gen,sel] = redundances[gen,sel] / gen.nevs
            }
          }
          red <- red + redundances[gen,sel]
        }
        
        # Score = relevance - relevance * redundance/ num of selected genes
        score <- rel - red/length(selected.genes) * rel
        
        # If actual score is max keep actual gene
        if (score > max.score){
          max.score <- score
          act.selected <- gen
        }
      }
      # Save gene with maximun score
      selected.genes[[act.selected]] <- max.score
      colnames(redundances)[length(selected.genes)] <- act.selected
    }
    return(selected.genes)
  }else{
    stop("The mode is unrecognized, please use mrmr or rf.")
  }
}


