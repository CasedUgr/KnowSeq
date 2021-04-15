#' DEGsExtraction performs the analysis to extract the Differentially Expressed Genes (DEGs) among the classes to compare.
#'
#' The function performs the analysis to extract the Differentially Expressed Genes (DEGs) among the classes to compare. The number of final DEGs can change depending on the p-value and the LFC indicated by parameters of the function. Furthermore, the function detects if the number of classes are greater than 2 to perform a multiclass DEGs analysis.
#' @param expressionMatrix The expressionMatrix parameter is an expression matrix or data.frame that contains the genes in the rows and the samples in the columns.
#' @param labels A vector or factors that contains the labels for each of the samples in the expressionMatrix parameter.
#' @param pvalue The value of the p-value which determines the DEGs. If one or more genes have a p-value lower or equal to the selected p-value, they would be considered as DEGs.
#' @param lfc The value of the LFC which determines the DEGs. If one or more genes have a LFC greater or equal to the selected LFC, they would be considered as DEGs.
#' @param cov This value only works when there are more than two classes in the labels. This parameter establishes a minimum number of pair of classes combination in which exists differential expression to consider a genes as expressed genes.
#' @param nmax This value only works when there are more than two classes in the labels. NMAX indicates the maximum number of DEGs selected for each class pair comparison.
#' @param multiDegsMethod Select the multiclass extraction method for the process: cov or nmax
#' @param number The maximum number of desired genes as output of limma. As default, the function returns all the extracted DEGs with the selected parameters.
#' @param CV A boolean value that has to be setted to TRUE if the user would to run a Cross-Validation DEGs extraction process.
#' @param numFolds This parameter indicates the number of folds for the Cross-Validation process.
#' @return A list that contains two objects. The table with statistics of the different DEGs and a reduced expression matrix which contains the DEGs and the samples.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation, genesNames = TRUE)
#'
#' DEGsInformation <- DEGsExtraction(expressionMatrix, labels, lfc = 2.0,
#' pvalue = 0.01, number = Inf)
#'
#' topTable <- DEGsInformation$Table
#'
#' DEGsMatrix <- DEGsInformation$DEGsMatrix

DEGsExtraction <- function(expressionMatrix, labels, pvalue=0.05, lfc = 1.0, cov = 1, nmax = 1, multiDegsMethod = "cov", number = Inf, CV = FALSE, numFolds = 5){

      if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
      if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
      if(!is.numeric(pvalue)){stop("The class of pvalue parameter must be numeric.")}
      if(!is.numeric(lfc)){stop("The class of lfc parameter must be numeric.")}
      if(!is.numeric(cov)){stop("The class of cov parameter must be numeric.")}
      if(!is.numeric(number)){stop("The class of number parameter must be numeric.")}

      if(is.character(labels)){ labels <- as.factor(labels) }
    
      cvDatasets <- list()
      cvDEGsResults <- list()
      cvDEGsList <- list()
        
      if (CV){
        
        cat("Applying DEGs extraction CV.\n")
        
        lengthValFold <- dim(expressionMatrix)[2]/numFolds
          
        positions <- rep(seq_len(dim(expressionMatrix)[2]))
        randomPositions <- sample(positions)
        expressionMatrix <- expressionMatrix[,randomPositions]
        labels <- labels[randomPositions]
        
        for(i in seq_len(numFolds)){
          
          valFold <- seq(round((i-1)*lengthValFold + 1 ), round(i*lengthValFold))
          foldData<- expressionMatrix[,valFold]
          foldLabels <- labels[valFold]
          cvDatasets[[i]] <- list(foldData,foldLabels)
      
        }

      } else{
          cvDatasets[[1]] <- list(expressionMatrix,labels)
          numFolds = 1
      }
  
      if(length(levels(labels)) == 2){

        cat("Two classes detected, applying limma biclass\n")
          
        for(i in seq_len(numFolds)){
          
          condition <- cvDatasets[[i]][[2]]
          DE.design <- model.matrix(~condition)
          fit <- lmFit(cvDatasets[[i]][[1]],DE.design)
          fit <- eBayes(fit)
          table <- topTable(fit, number = number, coef = 2, sort.by = "logFC", p.value = pvalue, adjust = "fdr", lfc = lfc)
          DEGsMatrix <- cvDatasets[[i]][[1]][rownames(table),]
          DEGsMatrix <- DEGsMatrix[unique(rownames(DEGsMatrix)),]
          
          cvDEGsResults[[i]] <- list(table,DEGsMatrix)
          cvDEGsList[[i]] <- rownames(table)
          
          names(cvDEGsResults[[i]]) <- c("DEGs_Table","DEGs_Matrix")
        }

        names(cvDEGsResults) <- rep(paste("Fold_",seq_len(numFolds),sep = ""),1)
        
        if (CV){
          commonDEGs <- Reduce(intersect,cvDEGsList)
          cat(paste("A total of ",length(commonDEGs)," common genes has been found\n", sep = ""))
          results <- list(cvDEGsResults, commonDEGs, expressionMatrix[commonDEGs,])
          names(results) <- c("CV_Results", "Common_DEGs", "Exprs_Matrix")
        }else{
          results <- list(cvDEGsResults[[1]])
          names(results) <- c("DEG_Results")  
        }
        

      }else if(length(levels(labels)) > 2){

        cat("More than two classes detected, applying limma multiclass\n")
        
        for(z in seq_len(numFolds)){
          
            condition <- cvDatasets[[z]][[2]]
            condition <- as.factor(condition)
            designMulti <- model.matrix(~0+condition)
            colnames(designMulti) = as.character(levels(condition))
            fitmicroMulti <- lmFit(cvDatasets[[z]][[1]], designMulti)
    
            contrasts <- as.character()
    
            for(i in c(seq_len(length(levels(condition))))){
              for(j in c(i+seq_len(length(levels(condition))))){
                if(j <= length(levels(condition))){
                  newContrast <- paste(levels(condition)[i],"-",levels(condition)[j],sep = "")
                  contrasts <- c(contrasts,newContrast)
                }
              }
            }
            
            if(CV){
              cat(paste("FOLD ",z,": \n",sep = ""))
            }
            
            cat(paste("Contrasts: ", contrasts,"\n",sep = ""))
            cat(table(condition))
            cat("\n")
            
            cont.matrixMulti = makeContrasts(contrasts = contrasts, levels=levels(condition))
            fitmicroMultiContrast = contrasts.fit(fitmicroMulti,cont.matrixMulti)
            fitmicroMultiContrast <- eBayes(fitmicroMultiContrast)
            
            if(multiDegsMethod == "cov"){
            res = decideTests(fitmicroMultiContrast,p.value=pvalue,lfc=lfc)
            ind = which(apply(res,1,function(x) {length(which(x != 0))>cov}) == TRUE)
    
              if(length(ind) > 0){
      
                lfcIndmatrix <- fitmicroMultiContrast$coefficients[ind,]
                multIndMatches <- res[ind,]
                multIndMatches <- abs(multIndMatches)
      
                DEGsMultiClass <-  cvDatasets[[z]][[1]][names(ind),]
                cvDEGsResults[[z]] <- list(fitmicroMultiContrast,lfcIndmatrix,DEGsMultiClass)
                cvDEGsList[[z]] <- rownames(DEGsMultiClass)
                names(cvDEGsResults[[z]]) <- c("DEGs_Table","MulticlassLFC","DEGs_Matrix")
                
      
              }else{
                stop("There are not genes that complains these restrictions, please change the p-value, lfc or cov.")
              }
            }else if(multiDegsMethod == "nmax"){
              
              res.val <- decideTests(fitmicroMultiContrast,p.value=pvalue,lfc = lfc)
              ind.val <- which(apply(res.val,1,function(x) {length(which(x != 0))>0}) == TRUE)
              lfcIndmatrix.sig <- fitmicroMultiContrast$coefficients[ind.val,]
              lfcIndmatrix.abs <- abs(lfcIndmatrix.sig)
              
              lfcs.sig <- as.data.frame(lfcIndmatrix.sig)
              lfcs.sig <- cbind(rownames(lfcs.sig),lfcs.sig)
              
              
              if(dim(lfcIndmatrix.abs)[1] >= nmax){
                genesSeveralMaxLFC <- rownames(lfcIndmatrix.abs)[apply(lfcIndmatrix.abs,2, order, decreasing = TRUE)]
                genesSeveralMaxLFC <- data.frame(matrix(unlist(genesSeveralMaxLFC), nrow = dim(lfcIndmatrix.abs)[1], byrow = FALSE), stringsAsFactors = FALSE)
                genesSeveralMaxLFC <- genesSeveralMaxLFC[1:nmax,]
                
                genesFilteredByLFC <- data.frame(matrix(character(),dim(genesSeveralMaxLFC)[1],dim(genesSeveralMaxLFC)[2]),stringsAsFactors = FALSE)
                
                for (i in 1:dim(genesSeveralMaxLFC)[1]){
                  for (j in 1:dim(genesSeveralMaxLFC)[2]){
                    if(lfcIndmatrix.abs[genesSeveralMaxLFC[i,j],j] >= lfc){
                      genesFilteredByLFC[i,j] <- as.character(genesSeveralMaxLFC[i,j])
                    }
                  }
                }
            
                genesFilteredByLFC <- as.data.frame(genesFilteredByLFC)
                colnames(genesFilteredByLFC) <- colnames(lfcIndmatrix.abs)
                genesSeveralMaxLFC <- unique(unlist(genesFilteredByLFC))
                genesSeveralMaxLFC <- genesSeveralMaxLFC[!is.na(genesSeveralMaxLFC)]
                
                DEGsMultiClass <-  cvDatasets[[z]][[1]][genesSeveralMaxLFC,]
                cvDEGsResults[[z]] <- list(fitmicroMultiContrast,lfcIndmatrix.abs,DEGsMultiClass)
                cvDEGsList[[z]] <- rownames(DEGsMultiClass)
                names(cvDEGsResults[[z]]) <- c("DEGs_Table","MulticlassLFC","DEGs_Matrix")
              }
            
            }else{
              stop("There are not genes that complains these restrictions, please change the p-value, lfc or nmax")
            }
            
        }
        
        names(cvDEGsResults) <- rep(paste("Fold_",seq_len(numFolds),sep = ""),1)
        
        if (CV){
          commonDEGs <- Reduce(intersect,cvDEGsList)
          cat(paste("A total of ",length(commonDEGs)," common genes has been found\n", sep = ""))
          results <- list(cvDEGsResults, commonDEGs, expressionMatrix[commonDEGs,])
          names(results) <- c("CV_Results", "Common_DEGs", "Exprs_Matrix")
        }else{
          results <- list(cvDEGsResults[[1]])
          names(results) <- c("DEG_Results")  
        }
        
      }

      invisible(results)

}
