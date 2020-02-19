#' getAnnotationFromEnsembl returns the required information about a list of genes from Ensembl biomart.
#'
#' The function returns the required information about a list of genes from Ensembl biomart. This list of genes can be Ensembl ID, gene names or either of the possible values admited by Ensembl biomart. Furthermore, the reference genome can be chosen depending on the necessity of the user.
#' @param values A list of genes that contains the names or IDs.
#' @param attributes A vector which contains the different information attributes that the Ensembl biomart admit.
#' @param filter The attribute used as filter to return the rest of the attributes.
#' @param notHSapiens A boolean value that indicates if the user wants the human annotation or another annotation available in BiomaRt. The possible not human dataset can be consulted by calling the following function: biomaRt::listDatasets(useMart("ensembl")).
#' @param notHumandataset A dataset identification from biomaRt::listDatasets(useMart("ensembl")).
#' @return A matrix that contains all the information asked to the attributes parameter.
#' @examples
#'myAnnotation <- getAnnotationFromEnsembl(c("KRT19","BRCA1"),attributes=c("ensembl_gene_id","percentage_gene_gc_content"),filter='external_gene_name',notHSapiens=FALSE)
#' myAnnotation <- getAnnotationFromEnsembl(c("MGP_129S1SvImJ_G0038602", "MGP_129S1SvImJ_G0007718"),attributes=c("percentage_gene_gc_content","ensembl_gene_id"),filter='ensembl_gene_id',notHSapiens = TRUE, notHumandataset = 'mm129s1svimj_gene_ensembl')


getAnnotationFromEnsembl <- function(values,attributes=c("ensembl_gene_id","external_gene_name","percentage_gene_gc_content"), filter="ensembl_gene_id", notHSapiens = FALSE, notHumandataset = ""){

  if(typeof(attributes) != "character"){stop("The parameter attributes must be a character vector that contains the wanted ensembl attributes.")}
  if(typeof(values) != "character"){stop("The parameter values must be a character vector that contains the genes IDs.")}
  if(typeof(filter) != "character"){stop("The parameter filter must be a character that contains one attribute used as filter.")}
  if(!is.logical(notHSapiens)){stop("notHSapiens parameter can only takes the values TRUE or FALSE.")}
  
  dir <- system.file("extdata", package="KnowSeq")
  dir <- 'inst/extdata'
  if(!notHSapiens){

      cat("Getting annotation of the Homo Sapiens...\n")
      cat(paste("Using reference genome 38.\n"))

      myAnnotation <- read.csv(paste(dir,"GRCh38Annotation.csv",sep = "/"))
      
      if (filter %in% colnames(myAnnotation) && length(intersect(attributes,colnames(myAnnotation))) == length(attributes)){
        if ( all(values != 'allGenome')) myAnnotation <- myAnnotation[myAnnotation[[filter]] %in% values,]
        myAnnotation <- myAnnotation[,union(attributes,filter)]
        return(myAnnotation)
      }
      else{
        dataset.name <- 'hsapiens_gene_ensembl'
        filename <- paste(dataset.name,'.csv',sep='')
      }
  }else{
    if(length(notHumandataset)[1] == 0 || is.null(notHumandataset)){
      
      stop("The notHumandataset is empty! Please, provide a right notHumandataset")
      
    }

    dataset.name <- notHumandataset
    filename <- paste(notHumandataset,'.csv',sep='')
  }


  cat(paste("Downloading annotation ", dataset.name, "...\n"))

  # Create query
  query = paste('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
                  <Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                  <Dataset name = "',dataset.name,'" interface = "default" >',sep='')
  if ( all(values != 'allGenome')) 
    query <- paste(query,'<Filter name = "',filter,'" value = "',sep='')
                  
  for ( value in values) query <- paste(query,value,',',sep='')
  query <- str_sub(query, 1, nchar(query)-1)
  query <- paste(query,'"/>',sep='')
  
  for (attribute in attributes)
    query <- paste(query,'<Attribute name = "',attribute,'" />',sep='')
  if (! filter %in% attributes ) 
    query <- paste(query,'<Attribute name = "',filter,'" />',sep='')
  query <- paste(query,'</Dataset></Query>',sep='')

  # Download annotation file
  request <- paste('http://www.ensembl.org/biomart/martservice?query=',query,sep='')
  download.file(request,method='wget',destfile = filename)

  myAnnotation <- read.csv(filename,header=FALSE)

  if( grepl('ERROR',myAnnotation[1,1]) ){
    
    stop('Error in query, please check attributes and filter')
    
  }

  colnames(myAnnotation) <- union(attributes,filter)

  system(paste('rm',filename))

  return(myAnnotation)
}




