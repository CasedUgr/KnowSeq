#' DEGsToDiseases obtains the information about what diseases are related to the DEGs indicated by parameter.
#'
#' The function obtains the information about what diseases are related to the DEGs indicated by parameter. For that, the function makes use of the web platforms gene2Diseases and targetValidation.
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param size The number of diseases to retrieve from targetValidation
#' @param disease Query a specific disease instead of retrieving the whole list of related diseases.
#' @param getEvidences Boolean. If true, for each gene, a list of found evidences for each disease will be returned.
#' @return A list which contains the information about the diseases associated to each genes or to a set of genes. If getEvidences is TRUE, found evidences for each case will be returned too.
#' @examples
#' diseases <- DEGsToDiseases(c("KRT19","BRCA1"), getEvidences = FALSE)
 
DEGsToDiseases <- function(geneList, size = 10, disease = "", getEvidences = FALSE){

  if(length(geneList)[1] == 0 || is.null(geneList)){

    stop("The geneList is empty! Please, provide a right geneList.")

  }
  
  ensembl_list <- getGenesAnnotation(geneList, attributes = c("ensembl_gene_id","external_gene_name"), filter = "external_gene_name")
  query_string = "
    query target($ensemblId: String!, $size: Int!){
        target(ensemblId: $ensemblId){
          id
          approvedSymbol
          approvedName
          biotype
          associatedDiseases(page: {index: 0, size: $size}){
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
  
  count_query_string = "
    query target($ensemblId: String!){
        target(ensemblId: $ensemblId){
          id
          approvedSymbol
          approvedName
          biotype
          associatedDiseases{
            count
          }
        }
      }
    "
  base = "https://api.platform.opentargets.org/api/v4/graphql"
  
  categories <- c("literature","rna_expression","genetic_association","somatic_mutation","known_drug","animal_model","affected_pathway")
  
  results <- list()
  ensembl_list <- ensembl_list[match(geneList, ensembl_list$external_gene_name),]
  
  if(nchar(disease) == 0){  
    cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")

    for(i in seq_len(length(ensembl_list$external_gene_name))){
      
      variables <- list("ensemblId" = ensembl_list$ensembl_gene_id[i],"size" = size)
      
      # Construct POST request body object with query string and variables
      post_body <- list(query = query_string, variables = variables)
      
      # Perform POST request
      r <- POST(url=base, body=post_body, encode='json')
      
      # Print data to RStudio console
      response <- content(r)$data
      nDiseases <- length(response$target$associatedDiseases$rows)
      
      if(!is.na(ensembl_list$ensembl_gene_id[i])){
        
          currentGene = matrix(0, nrow = nDiseases, ncol = 9)
          colnames(currentGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")
          
          #if (getEvidences) currentEvidences = vector("list", 0)
          
          for (j in seq_len(nDiseases)){
            
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
  } else{
    
    cat(paste("Obtaining scores for ", disease,"...\n", sep = ""))
    diseaseGene = matrix(0, nrow = length(ensembl_list$external_gene_name), ncol = 9)
    colnames(diseaseGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")

    for(i in seq_len(length(ensembl_list$external_gene_name))){
      
      variables <- list("ensemblId" = ensembl_list$ensembl_gene_id[i])
      
      post_body <- list(query = count_query_string, variables = variables)
      
      # Perform POST request
      r <- POST(url=base, body=post_body, encode='json')
      
      # Print data to RStudio console
      response <- content(r)$data
      
      count <- response$target$associatedDiseases$count
      
      # Construct POST request body object with query string and variables
      variables <- list("ensemblId" = ensembl_list$ensembl_gene_id[i], "size" = count)
      
      post_body <- list(query = query_string, variables = variables)
      
      # Perform POST request
      r <- POST(url=base, body=post_body, encode='json')
      
      # Print data to RStudio console
      response <- content(r)$data
      nDiseases <- length(response$target$associatedDiseases$rows)
      
      if(!is.na(ensembl_list$ensembl_gene_id[i])){
        
        currentGene = matrix(0, nrow = nDiseases, ncol = 9)
        colnames(currentGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")

        #if (getEvidences) currentEvidences = vector("list", 0)
        
        for (j in seq_len(nDiseases)){
          
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
          
          disPos <- which(tolower(currentGene[,1]) == tolower(disease))
          
          if(length(disPos) == 0){
            diseaseGene[i,] <- rep(0,9)
          }else{
            diseaseGene[i,] <- currentGene[disPos,]
          }
          
        }
        
        
      }else{
        cat(paste("Removing ", unique(geneList)[i], " from the list, there is no exists associated diseases for this gene.\n",sep = ""))
      }
    }
    
    rownames(diseaseGene) <- ensembl_list$external_gene_name
    results[[1]] <- list('summary'=diseaseGene)
    names(results)[1] <- "DiseaseGene"
    
    cat("Disease scores acquired successfully!\n")
    
  }
  
  results <- results[lengths(results) != 0]
  invisible(results)
}
