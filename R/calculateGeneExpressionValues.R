#' Calculates the gene expression values by using a matrix of counts from RNA-seq.
#'
#' Calculates the gene expression values by using a matrix of counts from RNA-seq. Furthermore, the conversion from Ensembl IDs to genes names is performed by default, but can be changed with the parameter genesNames.
#'
#' @param countsMatrix The original counts matrix returned by \code{\link{countsToMatrix}} function or a matrix with the gene Ensembl ID in the rows and the samples in the columns that contains the count values.
#' @param annotation A matrix that contains the Ensembl IDs, the gene name and the percentage gene gc content for the genes available in the expression matrix. This annotation could be extracted from the function \code{\link{getAnnotationFromEnsembl}}.
#' @param genesNames A boolean variable which indicates if the rownames of the expression matrix are the genes Names (Symbols) or the ensembl IDs.
#' @param notHuman A boolean variable which indicates if the gene length file is the default gene length human file or another file indicated by parameter.
#' @param notHumanGeneLengthCSV Path to the CSV file that contains the gene length of the specie to use.
#' @return A matrix that contains the gene expression values. The rownames are the genes names or the Ensembl IDs and the colnames are the samples.
#' @examples
#' dir <- system.file("extdata", package="KnowSeq")
#' load(paste(dir,"/expressionExample.RData",sep = ""))
#'
#' expressionMatrix <- calculateGeneExpressionValues(countsMatrix,myAnnotation, genesNames = TRUE)

calculateGeneExpressionValues <- function(countsMatrix,annotation,genesNames=TRUE,notHuman = FALSE, notHumanGeneLengthCSV = ""){

  if(!is.logical(genesNames)){stop("genesNames parameter can only takes the values TRUE or FALSE.")}
  if(!is.matrix(countsMatrix)){stop("The class of countsMatrix parameter must be matrix.")}
  if(!is.data.frame(annotation)){stop("The class of countsMatrix parameter must be data.frame.")}
  if(!is.logical(notHuman)){stop("notHuman parameter can only takes the values TRUE or FALSE.")}

  if(!notHuman){

    dir <- system.file("extdata", package="KnowSeq")
    geneLength <- read.csv(file.path(dir, "Genes_length_Homo_Sapiens.csv"), header = TRUE, sep = ",")
    geneLength <- geneLength[,c(-2,-3)]

  }else{

    if(file.exists(notHumanGeneLengthCSV)){cat("Not Human gene length CSV file found!\n")}else{stop("Not Human gene length CSV file not found, please revise the path to the file.\n")}
    geneLength <- read.csv(notHumanGeneLengthCSV, header = TRUE, sep = ",")

  }

  cat("Calculating gene expression values...\n")

  myGCannot <- annotation$percentage_gene_gc_content
  names(myGCannot) <- annotation$ensembl_gene_id
  myGCannot <- myGCannot[rownames(countsMatrix)]
  summary(myGCannot)

  mygenes <- intersect(rownames(countsMatrix),annotation$ensembl_gene_id)

  mylength <- setNames(geneLength[match(mygenes,geneLength$Gene_stable_ID), 2], nm = mygenes)
  
  if(any(is.na(mylength))){
    NaPos <- which(is.na(mylength) == TRUE)
    mylength <- mylength[-NaPos]

    mygenes <- mygenes[-NaPos]
    rownames <- annotation$external_gene_name[match(annotation$ensembl_gene_id[-NaPos], mygenes)]
  }else{
    
    rownames <- annotation$external_gene_name[match(annotation$ensembl_gene_id, mygenes)]
    
  }
  myGCannot <- myGCannot[match(mygenes,names(myGCannot))]

  mycqn <- cqn(countsMatrix[mygenes,], lengths = mylength, x = myGCannot, sizeFactors = apply(countsMatrix, 2, sum), verbose = TRUE)

  cqnValues <- mycqn$y + mycqn$offset
  expressionMatrix <- cqnValues - min(cqnValues) + 1

  if(genesNames){
    rownames(expressionMatrix) <- rownames
    expressionMatrix <- expressionMatrix[unique(rownames(expressionMatrix)),]
  }

  return(expressionMatrix)

}

