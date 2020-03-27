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
  
    
  cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
  base = "https://api.opentargets.io/v3/platform/public/association/filter?target="
  
  genePeti = character()
  results <- vector("list", 0)
  
  for(j in seq_len(length(unique(geneList)))){
    r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",geneList[j],"&size=1&filter=target",sep = ""))
    respon <- content(r_Ensembl)
    ensembl_id <- respon$data[[1]]$data$ensembl_gene_id
    
    if(!is.na(ensembl_id)){
      
      genePeti = paste(base, ensembl_id,"&direct=true&fields=association_score&fields=disease.efo_info.label&size=",size,sep = "")
      r <- GET(genePeti)
      response = content(r)
      if(response$size != 0){
          currentGene = matrix(, nrow = response$size, ncol = 9)
          if (getEvidences) currentEvidences = vector("list", 0)
          
          for(i in seq(response$size)){
            currentGene[i,1] = response$data[[i]]$disease$efo_info$label
            currentGene[i,2] = response$data[[i]]$association_score$overall
            currentGene[i,3] = response$data[[i]]$association_score$datatypes$literature
            currentGene[i,4] = response$data[[i]]$association_score$datatypes$rna_expression
            currentGene[i,5] = response$data[[i]]$association_score$datatypes$genetic_association
            currentGene[i,6] = response$data[[i]]$association_score$datatypes$somatic_mutation
            currentGene[i,7] = response$data[[i]]$association_score$datatypes$known_drug
            currentGene[i,8] = response$data[[i]]$association_score$datatypes$animal_model
            currentGene[i,9] = response$data[[i]]$association_score$datatypes$affected_pathway
  
            if (getEvidences) currentEvidences[[currentGene[i,1]]] <- DEGsEvidences(unique(geneList)[j],response$data[[i]]$disease$efo_info$label,verbose=FALSE)[[unique(geneList)[j]]]
          }
          colnames(currentGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")
          
          
          if (getEvidences) results[[j]] <- list('summary'=currentGene,'evidences'=currentEvidences)
          else results[[j]] <- list('summary'=currentGene)
          names(results)[j] <- unique(geneList)[j]
          
      }else{
        
        cat(paste("Removing ", unique(geneList)[j], " from the list, there is no exists associated diseases for this gene.\n",sep = ""))
        
      }
    }
  }
  
  cat("Diseases acquired successfully!\n")
  results <- results[lengths(results) != 0]
  invisible(results)
}
