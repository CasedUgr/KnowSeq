#' geneOntologyEnrichment obtains the information about what Gene Ontology terms are related to the DEGs.
#'
#' The function obtains the information about GO terms from the three differents ontologies that are related to the DEGs. The function also returns the description about each GO and a list of genes that are inside of each GO.
#' @param geneMatrix A matrix that contains the expression of the DEGs for each samples.
#' @param labels A vector that contains the labels of the samples of the DEGsMatrix.
#' @param identificator The identification methods for the genes. By default the identificator is the gene symbol or gene name.
#' @param mapping The annotation database to map the gene and GOs. By default is prepared for homo sapiens but can be changed for other species.
#' @param nGOs Maximun number of GOs to return of the total amount of top GOs. By default is equal to 10.
#' @param pvalCutOff The maximum p-value to considers that a genes is related with a GO term.
#' @return A list that contains a matrix for each of the possible ontologies and a matrix with the GOs for the three ontologies together.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' labelsGo <- gsub("Control",0,labels)
#' labelsGo <- gsub("Tumor",1,labelsGo)
#'
#' GOsList <- geneOntologyEnrichment(DEGsMatrix,labelsGo,nGOs = 20,pvalCutOff = 0.001)

geneOntologyEnrichment <- function(geneMatrix, labels, identificator = "SYMBOL", mapping = "org.Hs.eg.db",nGOs = 10, pvalCutOff = 0.01){


  if(!is.matrix(geneMatrix)){stop("The class of geneMatrix parameter must be matrix.")}
  if(!is.character(labels) && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
  if(!is.numeric(nGOs)){stop("The class of nGOs parameter must be numeric.")}
  if(!is.numeric(pvalCutOff)){stop("The class of pvalCutOff parameter must be numeric.")}

  PValueDiffGenes <-  getPvalues(geneMatrix, classlabel = as.integer(labels), alternative = "greater")

  cat("Searching GOs from BP Ontology...")
  GOdata <- new( "topGOdata", ontology="BP", allGenes = PValueDiffGenes, geneSel = function (allScore){return(allScore < pvalCutOff)}, annot=annFUN.org, mapping=mapping, ID = identificator)

  if(nGOs == Inf){

    nGOs =  length(genesInTerm(GOdata))

  }

  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher", pvalCutOff = pvalCutOff)
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks", pvalCutOff = pvalCutOff)
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks", pvalCutOff = pvalCutOff)

  allResBP <- GenTable(GOdata, classicFisher = resultFisher,  classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher",topNodes = nGOs)

  allGO = genesInTerm(GOdata)

  tableGO.DEGs <- allGO[allResBP$GO.ID]

  tableGO.DEGs.matrix.BP <- matrix(0L, nrow = length(tableGO.DEGs), ncol = 11)
  colnames(tableGO.DEGs.matrix.BP) <- c(colnames(allResBP),"GO_Genes","Description")

  base <- "https://www.ebi.ac.uk/QuickGO/services"

  for(i in seq_len(length(tableGO.DEGs))){

    tableGO.DEGs.matrix.BP[i,1] <- allResBP[i,1]
    tableGO.DEGs.matrix.BP[i,2] <- allResBP[i,2]
    tableGO.DEGs.matrix.BP[i,3] <- allResBP[i,3]
    tableGO.DEGs.matrix.BP[i,4] <- allResBP[i,4]
    tableGO.DEGs.matrix.BP[i,5] <- allResBP[i,5]
    tableGO.DEGs.matrix.BP[i,6] <- allResBP[i,6]
    tableGO.DEGs.matrix.BP[i,7] <- allResBP[i,7]
    tableGO.DEGs.matrix.BP[i,8] <- allResBP[i,8]
    tableGO.DEGs.matrix.BP[i,9] <- allResBP[i,9]
    petition <- paste("/ontology/go/terms/",gsub(":","%3A",names(tableGO.DEGs)[i]),sep = "")
    apiCall <- paste(base, petition,sep = "")
    get_GO <- httr::GET(apiCall)
    get_GO_text <- httr::content(get_GO, "text")
    get_GO_json <- fromJSON(get_GO_text, flatten = TRUE)
    tableGO.DEGs.matrix.BP[i,10] <- paste(as.character(tableGO.DEGs[[i]]),collapse=",")
    tableGO.DEGs.matrix.BP[i,11] <- as.character(get_GO_json$results$definition.text)

  }

  cat("Searching GOs from MF Ontology...")
  GOdata <- new( "topGOdata", ontology="MF", allGenes = PValueDiffGenes, geneSel = function (allScore){return(allScore < pvalCutOff)}, annot=annFUN.org, mapping=mapping, ID = identificator)

  if(nGOs == Inf){

    nGOs =  length(genesInTerm(GOdata))

  }

  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher", pvalCutOff = pvalCutOff)
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks", pvalCutOff = pvalCutOff)
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks", pvalCutOff = pvalCutOff)

  allResMF <- GenTable(GOdata, classicFisher = resultFisher,  classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher",topNodes = nGOs)

  allGO = genesInTerm(GOdata)

  tableGO.DEGs <- allGO[allResMF$GO.ID]

  tableGO.DEGs.matrix.MF <- matrix(0L, nrow = length(tableGO.DEGs), ncol = 11)
  colnames(tableGO.DEGs.matrix.MF) <- c(colnames(allResMF),"GO_Genes","Description")

  base <- "https://www.ebi.ac.uk/QuickGO/services"

  for(i in seq_len(length(tableGO.DEGs))){

    tableGO.DEGs.matrix.MF[i,1] <- allResMF[i,1]
    tableGO.DEGs.matrix.MF[i,2] <- allResMF[i,2]
    tableGO.DEGs.matrix.MF[i,3] <- allResMF[i,3]
    tableGO.DEGs.matrix.MF[i,4] <- allResMF[i,4]
    tableGO.DEGs.matrix.MF[i,5] <- allResMF[i,5]
    tableGO.DEGs.matrix.MF[i,6] <- allResMF[i,6]
    tableGO.DEGs.matrix.MF[i,7] <- allResMF[i,7]
    tableGO.DEGs.matrix.MF[i,8] <- allResMF[i,8]
    tableGO.DEGs.matrix.MF[i,9] <- allResMF[i,9]
    petition <- paste("/ontology/go/terms/",gsub(":","%3A",names(tableGO.DEGs)[i]),sep = "")
    apiCall <- paste(base, petition,sep = "")
    get_GO <- httr::GET(apiCall)
    get_GO_text <- httr::content(get_GO, "text")
    get_GO_json <- fromJSON(get_GO_text, flatten = TRUE)
    tableGO.DEGs.matrix.MF[i,10] <- paste(as.character(tableGO.DEGs[[i]]),collapse=",")
    tableGO.DEGs.matrix.MF[i,11] <- as.character(get_GO_json$results$definition.text)

  }

  cat("Searching GOs from CC Ontology...")
  GOdata <- new( "topGOdata", ontology="CC", allGenes = PValueDiffGenes, geneSel = function (allScore){return(allScore < pvalCutOff)}, annot=annFUN.org, mapping=mapping, ID = identificator)

  if(nGOs == Inf){

    nGOs =  length(genesInTerm(GOdata))

  }

  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher", pvalCutOff = pvalCutOff)
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks", pvalCutOff = pvalCutOff)
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks", pvalCutOff = pvalCutOff)

  allResCC <- GenTable(GOdata, classicFisher = resultFisher,  classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher",topNodes = nGOs)

  allGO = genesInTerm(GOdata)

  tableGO.DEGs <- allGO[allResCC$GO.ID]

  tableGO.DEGs.matrix.CC <- matrix(0L, nrow = length(tableGO.DEGs), ncol = 11)
  colnames(tableGO.DEGs.matrix.CC) <- c(colnames(allResCC),"GO_Genes","Description")

  base <- "https://www.ebi.ac.uk/QuickGO/services"

  for(i in seq_len(length(tableGO.DEGs))){

    tableGO.DEGs.matrix.CC[i,1] <- allResCC[i,1]
    tableGO.DEGs.matrix.CC[i,2] <- allResCC[i,2]
    tableGO.DEGs.matrix.CC[i,3] <- allResCC[i,3]
    tableGO.DEGs.matrix.CC[i,4] <- allResCC[i,4]
    tableGO.DEGs.matrix.CC[i,5] <- allResCC[i,5]
    tableGO.DEGs.matrix.CC[i,6] <- allResCC[i,6]
    tableGO.DEGs.matrix.CC[i,7] <- allResCC[i,7]
    tableGO.DEGs.matrix.CC[i,8] <- allResCC[i,8]
    tableGO.DEGs.matrix.CC[i,9] <- allResCC[i,9]
    petition <- paste("/ontology/go/terms/",gsub(":","%3A",names(tableGO.DEGs)[i]),sep = "")
    apiCall <- paste(base, petition,sep = "")
    get_GO <- httr::GET(apiCall)
    get_GO_text <- httr::content(get_GO, "text")
    get_GO_json <- fromJSON(get_GO_text, flatten = TRUE)
    tableGO.DEGs.matrix.CC[i,10] <- paste(as.character(tableGO.DEGs[[i]]),collapse=",")
    tableGO.DEGs.matrix.CC[i,11] <- as.character(get_GO_json$results$definition.text)

  }

  cat("Gene Ontology Enrichment done successfully!\n")

  allOntologiesTable <- rbind(tableGO.DEGs.matrix.BP,tableGO.DEGs.matrix.MF,tableGO.DEGs.matrix.CC)
  results <- list(tableGO.DEGs.matrix.BP,tableGO.DEGs.matrix.MF,tableGO.DEGs.matrix.CC,allOntologiesTable)
  names(results) <- c("BP Ontology GOs","MF Ontology GOs","CC Ontology GOs", "All Ontologies GO")
  invisible(results)

}
