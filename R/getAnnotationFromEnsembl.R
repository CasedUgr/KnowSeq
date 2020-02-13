#' getAnnotationFromEnsembl returns the required information about a list of genes from Ensembl biomart.
#'
#' The function returns the required information about a list of genes from Ensembl biomart. This list of genes can be Ensembl ID, gene names or either of the possible values admited by Ensembl biomart. Furthermore, the reference genome can be chosen depending on the necessity of the user.
#' @param values A list of genes that contains the names or IDs.
#' @param attributes A vector which contains the different information attributes that the Ensembl biomart admit.
#' @param filters The attributes used as filter to return the rest of the attributes.
#' @param notHSapiens A boolean value that indicates if the user wants the human annotation or another annotation available in BiomaRt. The possible not human dataset can be consulted by calling the following function: biomaRt::listDatasets(useMart("ensembl")).
#' @param notHumandataset A dataset identification from biomaRt::listDatasets(useMart("ensembl")).
#' @return A matrix that contains all the information asked to the attributes parameter.
#' @examples
#' myAnnotation <- getAnnotationFromEnsembl(c("ENSG00000210049","ENSG00000211459","ENSG00000210077"),notHSapiens=FALSE)
#' myAnnotation <- getAnnotationFromEnsembl(c("MGP_129S1SvImJ_G0038602", "MGP_129S1SvImJ_G0007718","MGP_129S1SvImJ_G0023218"),notHSapiens = TRUE, notHumandataset = 'mm129s1svimj_gene_ensembl')



getAnnotationFromEnsembl <- function(values,attributes=c("ensembl_gene_id","external_gene_name","percentage_gene_gc_content"), filters="ensembl_gene_id", notHSapiens = FALSE, notHumandataset = ""){

  if(typeof(attributes) != "character"){stop("The parameter attributes must be a character vector that contains the wanted ensembl attributes.")}
  if(typeof(values) != "character"){stop("The parameter values must be a character vector that contains the genes IDs.")}
  if(typeof(filters) != "character"){stop("The parameter filters must be a character vector that contains at least one attributes used as filter.")}
  if(!is.logical(notHSapiens)){stop("notHSapiens parameter can only takes the values TRUE or FALSE.")}

  if(!notHSapiens){

      cat("Getting annotation of the Homo Sapiens...\n")
      cat(paste("Using reference genome 38.\n"))

      myAnnotation <- read.csv('inst/extdata/GRCh38Annotation.csv')
      myAnnotation <- myAnnotation[myAnnotation$ensembl_gene_id %in% values,]
  }else{

    if(length(notHumandataset)[1] == 0 || is.null(notHumandataset)){stop("The notHumandataset is empty! Please, provide a right notHumandataset")}
    
    filename <- paste('inst/extdata/',notHumandataset,'.csv',sep='')

    if (file.exists(filename)){
      cat(paste("Getting annotation ", notHumandataset, "...\n"))
      myAnnotation <- read.csv(filename,header=TRUE)
    }
    else{
      cat(paste("Downloading annotation ", notHumandataset, "...\n"))
      query = paste('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
          <Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
          <Dataset name = "',notHumandataset,'" interface = "default" >',sep='')
      for (attribute in attributes)
        query <- paste(query,'<Attribute name = "',attribute,'" />',sep='')
      query <- paste(query,'</Dataset></Query>',sep='')
    
      request <- paste('http://www.ensembl.org/biomart/martservice?query=',query,sep='')
      download.file(request,method='wget',destfile = filename)
      myAnnotation <- read.csv(filename)
      
      # Set columns names as attributes
      colnames(myAnnotation) <- attributes
      write.csv(myAnnotation,filename,row.names=FALSE)

    }
    myAnnotation <- myAnnotation[myAnnotation$ensembl_gene_id %in% values,]
  }
  return(myAnnotation)
}



