#' DEGsEvidences function returns for each DEG a list of evidences that correlate it with the studied disease.
#'
#' DEGsEvidences function returns for each DEG a list of evidences that correlate it with the studied disease.
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param disease The name of a disease in order to obtain related evidences from target validation by using the DEGs indicated in the geneList parameter.
#' @param size The number of diseases to retrieve from targetValidation
#' @param verbose Boolean that indicates if progress messages are printed to stdout
#' @return A list which names are genes from geneList and which contains related evidences for each gene in geneList and indicated disease.
#' @examples
#' evidences <- DEGsEvidences(c("KRT19","BRCA1","TYMP"),'cancer')

DEGsEvidences <- function(geneList, disease,size=10, verbose=TRUE){
  if(length(geneList)[1] == 0 || is.null(geneList) ){

    stop("The geneList is empty! Please, provide a right geneList.")

  }
  if(disease == ''){
    
    stop("Please, indicate a disease name to acquire related evidences.")
    
  }
  
  ensembl_list <- getGenesAnnotation(geneList, attributes = c("ensembl_gene_id","external_gene_name"), filter = "external_gene_name")
  
  # Get disease id (it's necesary for evidences request)
  query_string_disease = "
  query target($disease: String!) {
	search(queryString: $disease){
    hits{
      id
      name
    }
  	total
	}
  }
  "
  base = "https://api.platform.opentargets.org/api/v4/graphql"

  variables <- list("disease" = disease)
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string_disease, variables = variables)
  
  # Perform POST request
  r <- POST(url=base, body=post_body, encode='json')
  
  # Print data to RStudio console
  response <- content(r)$data
    
  if (response$search$total == 0){
    stop("Disease not found")
  }
  
  disease.id <- response$search$hits[[1]]$id

  if (verbose) cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")

    # Create empty output
  info <- list()
  
  i=1
  # Iter in genes from geneList
  for(gene in unique(ensembl_list$ensembl_gene_id)){
    
    query_string_evidences = paste("
      query target($diseaseId: String=\"",disease.id,"\", $ensemblId: [String!]=\"",gene,"\",) {
    	disease(efoId: $diseaseId){
        name
        evidences(ensemblIds: $ensemblId){
          rows{
            literature
                  pathways{
                    id
                    name
                  }
                  drug{
                    id
                    name
                  }
          }
        }
      }
      }
    ",sep ="")
    
    # Construct POST request body object with query string and variables
    post_body <- list(query = query_string_evidences)
    
    # Perform POST request
    r <- POST(url=base, body=post_body, encode='json')
    
    # Print data to RStudio console
    response <- content(r)$data
    
    if( length(response) > 0 )ensembl_id <- gene
    else ensembl_id <- NULL
    
    if(!is.na(ensembl_id)){
      num_evidences <- length(response$disease$evidences$rows)
      currentGene = matrix(0, nrow = num_evidences, ncol = 3)
      colnames(currentGene) <- c("Literature", "Pathways", "Drugs")
      
      for(j in seq(num_evidences)){
        literature <- response$disease$evidences$rows[[j]]$literature[[1]]
        pathways <- response$disease$evidences$rows[[j]]$pathways[[1]]
        drug <- response$disease$evidences$rows[[j]]$drug
        if(!is.null(literature)){
          currentGene[[j,1]] <- literature
        }
        if(!is.null(pathways)){
          currentGene[[j,2]] <- pathways$name
        }
        if(!is.null(drug)){
          currentGene[[j,3]] <- drug$name
        }
      }
      info[[i]] <- list('summary'=currentGene)
      names(info)[i] <- unique(geneList)[i]
    }
    
    i <- i + 1
    
  }
  
  if (verbose) cat("Evidences acquired successfully!\n")
  invisible(info)
}