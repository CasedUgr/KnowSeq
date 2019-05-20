#' The function uses the DEGs to show graphically the expression of the samples in the pathways in which those genes appear.
#'
#' The function uses the DEGs to show graphically the expression of the samples in the pathways in which those genes appear. For that, the function makes use of a DEGsMatrix with the expression of the DEGs and the annotation of those DEGs in which appear the pathway or pathways of each DEGs. Internally, the function uses \code{\link{pathview}} to retrieve and colours the pathways, but a maximum number of 24 samples can be used. Furthermore, the function needs the expression matrix with all the genes in order to use them to colour the rest of the elements in the pathways.
#' @param DEGsMatrix A matrix that contains the expression of the DEGs for each samples. This matrix can be achieved by calling the function \code{\link{limmaDEGsExtraction}}. If the samples are more than 24, only the first 24 will be used to colour the pathways.
#' @param DEGsAnnotation A matrix that contains the gene names and the entrez IDs for the genes available in the DEGs matrix. This annotation can be obtained from the function \code{\link{getAnnotationFromEnsembl}}.
#' @param expressionMatrix A matrix that contains the expression of the all the genes available for each samples. If the samples are more than 24, only the first 24 will be used to colour the pathways.
#' @param expressionAnnotation A matrix that contains the gene names and the entrez IDs for all the genes available. This annotation can be obtained from the function \code{\link{getAnnotationFromEnsembl}}.
#' @param labels A vector that contains the labels of the samples for both the DEGsMatrix and the expressionMatrix.
#' @return Nothing to return.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' countsInfo <- read.csv(paste(dir,"/countFiles/mergedCountsInfo.csv",sep = ""))
#' 
#' countsInfo$Path <- paste(dir,"/countFiles/",countsInfo$Run,sep = "")
#' 
#' write.csv(countsInfo, file = "countsInfo.csv")
#'
#' countsInformation <- countsToMatrix("countsInfo.csv")
#'
#' countsMatrix <- countsInformation$countsMatrix
#' labels <- countsInformation$labels
#'
#' file.remove("countsInfo.csv")
#' 
#' myAnnotation <- getAnnotationFromEnsembl(rownames(countsMatrix),referenceGenome=37)
#'
#' expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation, genesNames = TRUE)
#'
#' DEGsInformation <- limmaDEGsExtraction(expressionMatrix, labels, lfc = 2.0,
#' pvalue = 0.01, number = Inf)
#'
#' topTable <- DEGsInformation$Table
#'
#' DEGsMatrix <- DEGsInformation$DEGsMatrix
#'
#' myDEGsAnnotation <- getAnnotationFromEnsembl(rownames(DEGsMatrix)[1:3],
#' referenceGenome=38,attributes = c("external_gene_name","entrezgene"),
#' filters = "external_gene_name")
#'
#' allMyAnnotation <-  getAnnotationFromEnsembl(rownames(expressionMatrix),
#'                                              referenceGenome=38,attributes = c("external_gene_name","entrezgene"),
#'                                              filters = "external_gene_name")
#'
#' DEGsPathwayVisualization(DEGsMatrix[1:3,], myDEGsAnnotation, expressionMatrix, allMyAnnotation, labels)

DEGsPathwayVisualization <- function(DEGsMatrix, DEGsAnnotation, expressionMatrix, expressionAnnotation, labels){

    if(!is.matrix(DEGsMatrix)){stop("The class of DEGsMatrix parameter must be matrix.")}
    if(!is.data.frame(DEGsAnnotation)){stop("The class of DEGsAnnotation parameter must be data.frame.")}
    if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
    if(!is.data.frame(expressionAnnotation)){stop("The class of expressionAnnotation parameter must be data.frame.")}

    DEGsAnnotation <- DEGsAnnotation[which(!is.na(DEGsAnnotation$entrezgene) == TRUE),]

    commonDEGs <- intersect(rownames(DEGsMatrix),unique(DEGsAnnotation$external_gene_name))
    posCommonDEGs <- match(rownames(DEGsMatrix[commonDEGs,]),DEGsAnnotation$external_gene_name)
    DEGsMatrix <- DEGsMatrix[commonDEGs,]
    rownames(DEGsMatrix) <- DEGsAnnotation$entrezgene[posCommonDEGs]

    expressionAnnotation <- expressionAnnotation[which(!is.na(expressionAnnotation$entrezgene) == TRUE),]

    commonGenes <- intersect(rownames(expressionMatrix),unique(expressionAnnotation$external_gene_name))
    posCommonGenes <- match(rownames(expressionMatrix[commonGenes,]),expressionAnnotation$external_gene_name)
    expressionMatrix <- expressionMatrix[commonGenes,]
    rownames(expressionMatrix) <- expressionAnnotation$entrezgene[posCommonGenes]

    pathways_unique <- character()

    cat("Retrieving DEGs associated pathways...\n")

    for(gene in DEGsAnnotation$entrezgene){

      get_GO <- httr::GET(paste("http://rest.kegg.jp/get/hsa:",gene,sep = ""))
      get_GO_text <- httr::content(get_GO, "text")
      pathways <- gsub('^.*PATHWAY\\s*|\\s*MODULE.*$', '', get_GO_text)
      pathways_unique <- c(pathways_unique,unique(as.character(unlist(str_extract_all(pathways,"hsa[a-zA-Z0-9]{5}")))))

    }

    pathways.unique <- unique(pathways_unique)
    pathways.unique <- pathways.unique[-which(pathways.unique == "hsa01100")]

    cat(paste("A total of ",length(pathways.unique)," pathways have been retrieved!\n", sep = ""))
    cat(paste(pathways.unique,"\n",sep = ""))
    pathways.dir <- paste(getwd(),"Pathways",sep = "/")
    dir.create(pathways.dir)

    expressionMatrixNorm <- expressionMatrix

    for(i in seq_len(dim(expressionMatrix)[1])){

      maxValue <- max(expressionMatrix[i,])
      minValue <- min(expressionMatrix[i,])
      expressionMatrixNorm[i,] <- ((expressionMatrix[i,] - minValue) / (maxValue - minValue)) * 2 - 1

    }

    if(dim(expressionMatrixNorm)[2] > 24){expressionMatrixNorm <- expressionMatrixNorm[,seq_len(24)]}

    for(pathway in pathways.unique){

      pathChecking <- httr::GET(paste("http://rest.kegg.jp/get/",pathway,sep = ""))

      if(pathChecking$status_code != 404){

        kegg.dir = paste(pathways.dir,pathway,sep = "/")
        dir.create(kegg.dir)

        tryCatch(
          {

            pathways.out <- pathview(gene.data = expressionMatrixNorm, pathway.id = pathway, species = "hsa", new.signature=FALSE)
            fileMove(from = paste(pathway,".png",sep = ""),
                     to =  paste(kegg.dir,"/",pathway,".png",sep = ""))
            fileMove(from = paste(pathway,".pathview.png",sep = ""),
                     to =  paste(kegg.dir,"/",pathway,".pathview.png",sep = ""))
            fileMove(from = paste(pathway,".pathview.multi.png",sep = ""),
                     to =  paste(kegg.dir,"/",pathway,".pathview.multi.png",sep = ""))
            fileMove(from = paste(pathway,".xml",sep = ""),
                     to =  paste(kegg.dir,"/",pathway,".xml",sep = ""))

          },
          error=function(cond) {
            cat("This pathway cannot be processed successfully...")
            file.remove(paste(pathway,"*"))
            }
          )


      }else{

        cat(paste("Pathway ", pathway, " cannot be downloaded due to a 404 error.\n",sep = ""))

      }
    }

}
