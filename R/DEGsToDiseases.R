#' DEGsToDiseases obtains the information about what diseases are related to the DEGs indicated by parameter.
#'
#' The function obtains the information about what diseases are related to the DEGs indicated by parameter. For that, the function makes use of the web platforms gene2Diseases and targetValidation.
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param minCitation Minimum number of citations of each genes in a disease to consider the genes related with the disease.
#' @param size The number of diseases to retrieve from targetValidation
#' @param getEvidences Boolean. If true, for each gene, a list of found evidences for each disease will be returned.
#' @return A list which contains the information about the diseases associated to each genes or to a set of genes. If getEvidences is TRUE, found evidences for each case will be returned too.
#' @examples
#' diseases <- DEGsToDiseases(c("KRT19","BRCA1"),getEvidences = FALSE)

DEGsToDiseases <- function(geneList, minCitation = 5, size = 10, getEvidences = FALSE){

  if(length(geneList)[1] == 0 || is.null(geneList)){

    stop("The geneList is empty! Please, provide a right geneList.")

  }
  
  ensembl_list <- getGenesAnnotation(geneList, attributes = c("ensembl_gene_id","external_gene_name"), filter = "external_gene_name")
    
  cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
  
  # Build query string
  query_string = "
  query target($ensemblId: String!){
      target(ensemblId: $ensemblId){
        id
        approvedSymbol
        approvedName
        bioType
        associatedDiseases{
          count
          rows{
            score
            datatypeScores{
              id
              score
            }
            disease{
              name
            }
          }
        }
      }
    }
  "
  
  base = "https://api.platform.opentargets.org/api/v4/graphql"
  
  categories <- c("literature","rna_expression","genetic_association","somatic_mutation","known_drug","animal_model","affected_pathway")
  
  results <- list()
  
  for(i in seq_len(length(unique(ensembl_list)))){
    
    variables <- list("ensemblId" = ensembl_list$ensembl_gene_id[i])
    
    # Construct POST request body object with query string and variables
    post_body <- list(query = query_string, variables = variables)
    
    # Perform POST request
    r <- POST(url=base, body=post_body, encode='json')
    
    # Print data to RStudio console
    response <- content(r)$data
    
    if(!is.na(ensembl_list$ensembl_gene_id[i])){
      
        currentGene = matrix(0, nrow = 25, ncol = 9)
        colnames(currentGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")
        
        #if (getEvidences) currentEvidences = vector("list", 0)
        
        for (j in seq_len(25)){
          
          currentGene[j,1] <- response$target$associatedDiseases$rows[[j]]$disease$name
          currentGene[j,2] <- response$target$associatedDiseases$rows[[j]]$score
          
          scores <- length(response$target$associatedDiseases$rows[[j]]$datatypeScores)
          
          for (k in seq_len(scores)){
            
            indexCategory <- match(response$target$associatedDiseases$rows[[j]]$datatypeScores[[k]]$id, categories) + 2
            currentGene[j,indexCategory] <- response$target$associatedDiseases$rows[[j]]$datatypeScores[[k]]$score
            
          }
          #if (getEvidences) currentEvidences[[currentGene[i,1]]] <- DEGsEvidences(unique(geneList)[j],response$data[[i]]$disease$efo_info$label,verbose=FALSE)[[unique(geneList)[j]]]
    
        #if (getEvidences) results[[j]] <- list('summary'=currentGene,'evidences'=currentEvidences)
        #else 
          results[[i]] <- list('summary'=currentGene)
          names(results)[i] <- unique(geneList)[i]
        }    
      }else{
        
        cat(paste("Removing ", unique(geneList)[i], " from the list, there is no exists associated diseases for this gene.\n",sep = ""))
        
      }
    }
  
  cat("Diseases acquired successfully!\n")
  results <- results[lengths(results) != 0]
  invisible(results)
}
