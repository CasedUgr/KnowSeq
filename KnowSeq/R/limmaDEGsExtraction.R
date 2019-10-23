#' limmaDEGsExtraction performs the analysis to extract the Differentially Expressed Genes (DEGs) among the classes to compare.
#'
#' The function performs the analysis to extract the Differentially Expressed Genes (DEGs) among the classes to compare. The number of final DEGs can change depending on the p-value and the LFC indicated by parameters of the function. Furthermore, the function detects if the number of classes are greater than 2 to perform a multiclass DEGs analysis.
#' @param expressionMatrix The expressionMatrix parameter is an expression matrix or data.frame that contains the genes in the rows and the samples in the columns.
#' @param labels A vector or factors that contains the labels for each of the samples in the expressionMatrix parameter.
#' @param pvalue The value of the p-value which determines the DEGs. If one or more genes have a p-value lower or equal to the selected p-value, they would be considered as DEGs.
#' @param lfc The value of the LFC which determines the DEGs. If one or more genes have a LFC greater or equal to the selected LFC, they would be considered as DEGs.
#' @param cov This value only works when there are more than two classes in the labels. This parameter stablishs a minimum number of pair of classes combination in which exists differential expression to consider a genes as expressed genes.
#' @param number The maximum number of desired genes as output of limma. As default, the function returns all the extracted DEGs with the selected parameters.
#' @param svaCorrection A logical variable that represents if the model for limma is calculated or indicated by parameter from the output of \code{\link{batchEffectRemoval}} function by using sva method.
#' @param svaMod The model calculated by \code{\link{batchEffectRemoval}} function by using sva method.
#' @return A list that contains two objects. The table with statistics of the different DEGs and a reduced expression matrix which contains the DEGs and the samples.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation, genesNames = TRUE)
#'
#' DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = 2.0,
#' pvalue = 0.01, number = Inf)
#'
#' topTable <- DEGsInformation$Table
#'
#' DEGsMatrix <- DEGsInformation$DEGsMatrix

limmaDEGsExtraction <- function(expressionMatrix, labels, pvalue=0.05, lfc = 1.0, cov = 1,number = Inf, svaCorrection = FALSE, svaMod){

      if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
      if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
      if(!is.logical(svaCorrection)){stop("svaCorrection parameter can only takes the values TRUE or FALSE.")}
      if(!is.numeric(pvalue)){stop("The class of pvalue parameter must be numeric.")}
      if(!is.numeric(lfc)){stop("The class of lfc parameter must be numeric.")}
      if(!is.numeric(cov)){stop("The class of cov parameter must be numeric.")}
      if(!is.numeric(number)){stop("The class of number parameter must be numeric.")}

      if(is.character(labels)){ labels <- as.factor(labels) }

      if(length(levels(labels)) == 2){

          cat("Two classes detected, applying limma biclass\n")
          if(!svaCorrection){

            condition <- labels
            DE.design <- model.matrix(~condition)
            fit <- lmFit(expressionMatrix,DE.design)

          }else{

            fit <- lmFit(expressionMatrix,svaMod)

          }

          fit <- eBayes(fit)
          table <- topTable(fit, number = number, coef = 2, sort.by = "logFC", p.value = pvalue, adjust = "fdr", lfc = lfc)
          DEGsMatrix <- expressionMatrix[rownames(expressionMatrix) %in% rownames(table),]
          DEGsMatrix <- DEGsMatrix[unique(rownames(DEGsMatrix)),]

          results <- list(table,DEGsMatrix)
          names(results) <- c("Table","DEGsMatrix")

      }else if(length(levels(labels)) > 2){

        cat("More than two classes detected, applying limma multiclass\n")

        condition <- labels
        designMulti <- model.matrix(~0+condition)
        colnames(designMulti) = levels(condition)
        fitmicroMulti <- lmFit(expressionMatrix, designMulti)

        contrasts <- as.character()

        for(i in c(seq_len(length(levels(condition))))){
          for(j in c(i+seq_len(length(levels(condition))))){
            if(j <= length(levels(condition))){
              newContrast <- paste(levels(condition)[i],"-",levels(condition)[j],sep = "")
              contrasts <- c(contrasts,newContrast)
            }
          }
        }

        cat(paste("Contrasts: ", contrasts,"\n",sep = ""))

        cont.matrixMulti = makeContrasts(contrasts = contrasts, levels=levels(condition))
        fitmicroMultiContrast = contrasts.fit(fitmicroMulti,cont.matrixMulti)
        fitmicroMultiContrast <- eBayes(fitmicroMultiContrast)
        res = decideTests(fitmicroMultiContrast,p.value=pvalue,lfc=lfc)
        ind = which(apply(res,1,function(x) {length(which(x != 0))>cov}) == TRUE)

        if(length(ind) > 0){

          lfcIndmatrix <- fitmicroMultiContrast$coefficients[ind,]
          multIndMatches <- res[ind,]
          multIndMatches <- abs(multIndMatches)

          DEGsMultiClass <- expressionMatrix[names(ind),]
          DEGsMultiClass <- t(DEGsMultiClass)
          DEGsMultiClass <- cbind(DEGsMultiClass,labels)

          results <- list(fitmicroMultiContrast,DEGsMultiClass)
          names(results) <- c("Table","DEGsMatrix")


        }else{
          stop("There are not genes that complains these restrictions, please change the p-value, lfc or cov.")
        }

      }

      invisible(results)

}
