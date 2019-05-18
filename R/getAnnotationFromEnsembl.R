#' getAnnotationFromEnsembl returns the required information about a list of genes from Ensembl biomart.
#'
#' The function returns the required information about a list of genes from Ensembl biomart. This list of genes can be Ensembl ID, gene names or either of the possible values admited by Ensembl biomart. Furthermore, the reference genome can be chosen depending on the necessity of the user.
#' @param values A list of genes that contains the names or IDs.
#' @param attributes A vector which contains the different information attributes that the Ensembl biomart admit.
#' @param filters The attributes used as filter to return the rest of the attributes.
#' @param referenceGenome The human reference genome used to return the annotation. The possibilities are two: 37 and 38
#' @param notHSapiens A boolean value that indicates if the user wants the human annotation or another annotation available in BiomaRt. The possible not human dataset can be consulted by calling the following function: biomaRt::listDatasets(useMart("ensembl")).
#' @param notHumandataset A dataset identification from biomaRt::listDatasets(useMart("ensembl")).
#' @return A matrix that contains all the information asked to the attributes parameter.
#' @examples
#' myAnnotation <- getAnnotationFromEnsembl(c("ENSG00000210049","ENSG00000211459","ENSG00000210077"),referenceGenome=37)

getAnnotationFromEnsembl <- function(values,attributes=c("ensembl_gene_id","external_gene_name","percentage_gene_gc_content","gene_biotype"), filters="ensembl_gene_id", referenceGenome = 38, notHSapiens = FALSE, notHumandataset = ""){

  if(typeof(attributes) != "character"){stop("The parameter attributes must be a character vector that contains the wanted ensembl attributes.")}
  if(typeof(values) != "character"){stop("The parameter values must be a character vector that contains the genes IDs.")}
  if(typeof(filters) != "character"){stop("The parameter filters must be a character vector that contains at least one attributes used as filter.")}
  if(!is.logical(notHSapiens)){stop("notHSapiens parameter can only takes the values TRUE or FALSE.")}

  if(!notHSapiens){

      cat("Downloading annotation of the Homo Sapiens...\n")

      if(referenceGenome == 37){
          cat(paste("Using reference genome 37.\n"))
          mart  <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
                                    path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
      }else if(referenceGenome == 38){
        cat(paste("Using reference genome 38.\n"))
        mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
      }else{
        stop("The version of the reference genome is wrong. Please use 37 or 38.")
      }

      myAnnotation <- getBM(attributes = attributes, filters = filters, values = values, mart = mart)

  }else{

    if(length(notHumandataset)[1] == 0 || is.null(notHumandataset)){stop("The notHumandataset is empty! Please, provide a right notHumandataset")}

    cat(paste("Downloading annotation ", notHumandataset, "...\n"))
    mart <- useMart("ensembl", dataset=notHumandataset)

    myAnnotation <- getBM(attributes = attributes, filters = filters, values = values, mart = mart)

  }

  return(myAnnotation)

}
