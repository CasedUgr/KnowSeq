
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param disease The name of a disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @param minCitation Minimum number of citations of each genes in a disease to consider the genes related with the disease.
#' @param size The number of diseases to retrieve from targetValidation
#' @return A list which names are genes from geneList and which contains related evidences for each gene in geneList and indicated disease.

DEGsEvidences <- function(geneList, disease, minCitation = 5, size = 10){
  if(length(geneList)[1] == 0 || is.null(geneList) ){
    
    stop("The geneList is empty! Please, provide a right geneList.")
    
  }
  if(disease == ''){
    
    stop("Please, indicate a disease name to acquire related evidences.")
    
  }
  
  # Get disease id
  r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease,"&size=1&filter=disease",sep = ""))
  respon <- content(r_Ensembl)
  if (respon$size == 0){
    stop("Disease not found")
  }
  disease.id <- respon$data[[1]]$id

  cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
  base = "https://api.opentargets.io/v3/platform/public/evidence/filter?target="
  
  genePeti = character()
  results <- vector("list", 0)
  info <- list()
  
  for(j in seq_len(length(unique(geneList)))){
    # Get gene id
    r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",geneList[j],"&size=1&filter=target",sep = ""))
    respon <- content(r_Ensembl)
    if( length(respon$data) > 0 )ensembl_id <- respon$data[[1]]$data$ensembl_gene_id
    else ensembl_id <- NULL
    
    
    if(length(ensembl_id) > 0){
      disease.req <- paste(base,ensembl_id,"&disease=",disease.id,"&direct=true&size=100",sep = "")
      disease.req <- GET(disease.req)
      response.disease = content(disease.req)
      
      if(response.disease$size != 0){
        act.evidences = matrix(, nrow = response.disease$size, ncol = 6)
        
        for(k in seq(response.disease$size)){
          if (response.disease$data[[k]]$target$gene_info$geneid == ensembl_id){
            type <- response.disease$data[[k]]$type
            act.evidences[k,1] = type
            act.evidences[k,2] = response.disease$data[[k]]$evidence$evidence_codes[[1]]
            act.evidences[k,3] = response.disease$data[[k]]$scores$association_score

            if (type == 'known_drug'){
              act.evidences[k,4] = response.disease$data[[k]]$drug$molecule_name
              act.evidences[k,5] = response.disease$data[[k]]$drug$molecule_type
            }
            else if (type == 'literature'){
              act.evidences[k,4] = response.disease$data[[k]]$evidence$literature_ref$lit_id # url
            }
            else if(type == 'affected_pathway'){
              if (length(response.disease$data[[k]]$evidence$resource_score$method$reference) > 0){
                act.evidences[k,4] = response.disease$data[[k]]$evidence$resource_score$method$reference # url
              }
              if (length(response.disease$data[[k]]$evidence$known_mutations)>0){
                act.evidences[k,5] = response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name
              }
              if(length( response.disease$data[[k]]$evidence$urls[[1]]$url ) > 0){
                act.evidences[k,6] = response.disease$data[[k]]$evidence$urls[[1]]$url # reactome url
              }
            }
            else if(type == 'animal_model'){
              act.evidences[k,4] = response.disease$data[[k]]$evidence$biological_model$is_associated
              act.evidences[k,5] = response.disease$data[[k]]$evidence$biological_model$species
            }
            else if (type == 'rna_expression'){
              if(length(response.disease$data[[k]]$literature$references[[1]]$lit_id>0))
                act.evidences[k,4] = response.disease$data[[k]]$literature$references[[1]]$lit_id[1]
              act.evidences[k,5] = response.disease$data[[k]]$evidence$comparison_name # ej: "'oral squamous cell carcinoma' vs 'normal'"
            }
            else if (type == 'genetic_association'){
              if(length(response.disease$data[[k]]$literature$references[[1]]$list_id > 0))
                act.evidences[k,4] = response.disease$data[[k]]$literature$references[[1]]$lit_id # paper
              if(length(response.disease$data[[k]]$evidence$gene2variant$functional_consequence)>0)
                act.evidences[k,5] = response.disease$data[[k]]$evidence$gene2variant$functional_consequence # url
            }
            else if (type == 'somatic_mutation'){
              act.evidences[k,4] = response.disease$data[[k]]$evidence$known_mutations[[1]]$functional_consequence #url
              act.evidences[k,5] = response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name # ej: "sequence_alteration"
              if (length(response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer) > 0)
                act.evidences[k,6] = response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer # ej: "oncogene, TSG"
            }
            else{
              print(paste('Unknown evidence type',response.disease$data[[k]]$type))
            }
          }
        }
        info[[geneList[j]]] <- act.evidences
      }else{
        info[[geneList[j]]] <- rbind(c('literature',rep('*',5)),c('literature',rep('*',5)))
        cat(paste( unique(geneList)[j], " isn't associated with ",disease,".\n",sep = ""))
      }
    }else{
      info[[geneList[j]]] <- rbind(c('literature',rep('*',5)),c('literature',rep('*',5)))
    }
  }
  
  cat("Evidences acquired successfully!\n")
  results <- results[lengths(results) != 0]
  invisible(info)
}