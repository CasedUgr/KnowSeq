#' The function uses the DEGs to show graphically the expression of the samples in the pathways in which those genes appear.
#'
#' The function uses the DEGs to show graphically the expression of the samples in the pathways in which those genes appear. For that, the function makes use of a DEGsMatrix with the expression of the DEGs and the annotation of those DEGs in which appear the pathway or pathways of each DEGs. Internally, the function uses \code{\link{pathview}} to retrieve and colours the pathways, but a maximum number of 24 samples can be used. Furthermore, the function needs the expression matrix with all the genes in order to use them to colour the rest of the elements in the pathways.
#' @param DEGsMatrix A matrix that contains the expression of the DEGs for each samples. This matrix can be achieved by calling the function \code{\link{DEGsExtraction}}. If the samples are more than 24, only the first 24 will be used to colour the pathways.
#' @param DEGsAnnotation A matrix that contains the gene names and the entrez IDs for the genes available in the DEGs matrix. This annotation can be obtained from the function \code{\link{getGenesAnnotation}}.
#' @param expressionMatrix A matrix that contains the expression of the all the genes available for each samples. If the samples are more than 24, only the first 24 will be used to colour the pathways.
#' @param expressionAnnotation A matrix that contains the gene names and the entrez IDs for all the genes available. This annotation can be obtained from the function \code{\link{getGenesAnnotation}}.
#' @param labels A vector that contains the labels of the samples for both the DEGsMatrix and the expressionMatrix.
#' @return Nothing to return.
#' @examples
#' \dontrun{DEGsPathwayVisualization(DEGsMatrix, myDEGsAnnotation, expressionMatrix, allMyAnnotation, labels)}
#' 

DEGsPathwayVisualization <- function(DEGsMatrix, DEGsAnnotation, expressionMatrix, expressionAnnotation, labels){

    if(!is.matrix(DEGsMatrix)){stop("The class of DEGsMatrix parameter must be matrix.")}
    if(!is.data.frame(DEGsAnnotation)){stop("The class of DEGsAnnotation parameter must be data.frame.")}
    if(!is.matrix(expressionMatrix)){stop("The class of expressionMatrix parameter must be matrix.")}
    if(!is.data.frame(expressionAnnotation)){stop("The class of expressionAnnotation parameter must be data.frame.")}

    DEGsAnnotation <- DEGsAnnotation[which(!is.na(DEGsAnnotation$entrezgene_id) == TRUE),]

    commonDEGs <- intersect(rownames(DEGsMatrix),unique(DEGsAnnotation$external_gene_name))
    posCommonDEGs <- match(rownames(DEGsMatrix[commonDEGs,]),DEGsAnnotation$external_gene_name)
    DEGsMatrix <- DEGsMatrix[commonDEGs,]
    rownames(DEGsMatrix) <- DEGsAnnotation$entrezgene_id[posCommonDEGs]

    expressionAnnotation <- expressionAnnotation[which(!is.na(expressionAnnotation$entrezgene_id) == TRUE),]

    commonGenes <- intersect(rownames(expressionMatrix),unique(expressionAnnotation$external_gene_name))
    posCommonGenes <- match(rownames(expressionMatrix[commonGenes,]),expressionAnnotation$external_gene_name)
    expressionMatrix <- expressionMatrix[commonGenes,]
    rownames(expressionMatrix) <- expressionAnnotation$entrezgene_id[posCommonGenes]

    pathways_unique <- character()

    cat("Retrieving DEGs associated pathways...\n")

    for(gene in DEGsAnnotation$entrezgene_id){

      if(!is.na(gene)){
          get_GO <- GET(paste("http://rest.kegg.jp/get/hsa:",gene,sep = ""))
          get_GO_text <- content(get_GO, "text")
          pathway_start <- str_locate_all(pattern = "PATHWAY", get_GO_text)[[1]][2]
          pathway_end <- str_locate_all(pattern = "BRITE", get_GO_text)[[1]][1]
          pathways <- substr(get_GO_text,pathway_start+1,pathway_end-1)
          pathways_unique <- c(pathways_unique,unique(as.character(unlist(str_extract_all(pathways,"hsa[a-zA-Z0-9]{5}")))))
      }
    }

    naPos <- which(is.na(pathways_unique) == TRUE)
    if(length(naPos) != 0)
      pathways_unique <- pathways_unique[-naPos]
    
    pathways.unique <- unique(pathways_unique)

    cat(paste("A total of ",length(pathways.unique)," pathways have been retrieved!\n", sep = ""))
    cat(paste(pathways.unique,"\n",sep = ""))
    pathways.dir <- paste(getwd(),"Pathways",sep = "/")
    dir.create(pathways.dir)

    expressionMatrixNorm <- expressionMatrix

    expressionMatrixNorm = vapply(as.data.frame(t(expressionMatrixNorm)), function(x){ 
      max = max(x)
      min = min(x)
      x = ((x-min)/(max-min))*2-1}, double(ncol(expressionMatrixNorm)), USE.NAMES = TRUE)
    
    expressionMatrixNorm <- t(expressionMatrixNorm)
    colnames(expressionMatrixNorm) <- colnames(expressionMatrix)
    
    if(dim(expressionMatrixNorm)[2] > 24){expressionMatrixNorm <- expressionMatrixNorm[,seq_len(24)]}

    for(pathway in pathways.unique){

      pathChecking <- GET(paste("http://rest.kegg.jp/get/",pathway,sep = ""))

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
