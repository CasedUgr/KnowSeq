#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param disease The name of a disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @param subdisease The name of a particular subdisease from disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter. Default '' 
#' @param minCitation Minimum number of citations of each genes in a disease to consider the genes related with the disease.
#' @param size The number of diseases to retrieve from targetValidation
#' @return A list which names are genes from geneList and which contains related evidences for each gene in geneList and indicated disease.
#' @example
#' evidences <- DEGsEvidences(c("KRT19","BRCA1","TYMP"),'cancer')

DEGsEvidences <- function(geneList, disease, subdisease='', minCitation = 5, size = 10){
  if(length(geneList)[1] == 0 || is.null(geneList) ){
    
    stop("The geneList is empty! Please, provide a right geneList.")
    
  }
  if(disease == ''){
    
    stop("Please, indicate a disease name to acquire related evidences.")
    
  }
  # Get disease id (it's necesary for evidences request)
  r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease,"&size=1&filter=disease",sep = ""))
  respon <- content(r_Ensembl)
  
  if ( 'size' %in% names(respon) && respon$size == 0){
    stop("Disease not found")
  }
  disease.id <- respon$data[[1]]$id
  
  cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
  base = "https://api.opentargets.io/v3/platform/public/evidence/filter?target="
  
  # Create empty output
  info <- list()
  
  # Iter in genes from geneList
  for(j in seq(length(unique(geneList)))){
    # Get gene id (it's necesary for evidences request)
    r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",geneList[j],"&size=1&filter=target",sep = ""))
    respon <- content(r_Ensembl)
    if( length(respon$data) > 0 )ensembl_id <- respon$data[[1]]$data$ensembl_gene_id
    else ensembl_id <- NULL
    
    # If gene is found
    if(length(ensembl_id) > 0){
      disease.req <- paste(base,ensembl_id,"&disease=",disease.id,"&direct=true&size=",size,sep = "")
      disease.req <- GET(disease.req)
      response.disease = content(disease.req)
      
      # If response is not empty
      if(response.disease$size != 0){
        # Empty matrix that will contain information from obtained evidences for actual gene
        evidences = list()
        
        # Iter in found evidences
        for(k in seq(response.disease$size)){
          # Check if disease is matching
          if ( (subdisease=='' && grepl(disease,response.disease$data[[k]]$disease$efo_info$label))
               || subdisease == response.disease$data[[k]]$disease$efo_info$label){
            # Create empty evidence
            act.evidence <- list()
            
            # Check if gene is matching
            if (response.disease$data[[k]]$target$gene_info$geneid == ensembl_id){
              type <- response.disease$data[[k]]$type
              # Save evidence score
              act.evidence <- list('score' = response.disease$data[[k]]$scores$association_score,
                                   'evidence' = c())
              # Save evidence codes
              act.evidence['code'] <- c()
              for (code.info in response.disease$data[[k]]$evidence$evidence_codes_info)
                act.evidence['codes'] <- c(code.info[[1]]$label)
              
              # Save information depending on evidence type
              if (type == 'known_drug'){
                act.evidence$evidence <- c(act.evidence$evidence,
                                           response.disease$data[[k]]$drug$molecule_name, # Drug name
                                           response.disease$data[[k]]$drug$molecule_type) # Molecule type
              }
              else if (type == 'literature'){
                act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$literature_ref$lit_id) # paper url
              }
              else if(type == 'affected_pathway'){
                if (length(response.disease$data[[k]]$evidence$resource_score$method$reference) > 0){
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$resource_score$method$reference) # paper url
                } else act.evidence$evidence <- c(act.evidence$evidence,'*')
                if (length(response.disease$data[[k]]$evidence$known_mutations)>0){
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name)
                } else act.evidence$evidence <- c(act.evidence$evidence,'*')
                if(length( response.disease$data[[k]]$evidence$urls[[1]]$url ) > 0){
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$urls[[1]]$url) # reactome url
                } else act.evidence$evidence <- c(act.evidence$evidence,'*')
              }
              else if(type == 'animal_model'){
                act.evidence$evidence <- c(act.evidence$evidence,
                                           response.disease$data[[k]]$evidence$biological_model$is_associated, # Boolean
                                           response.disease$data[[k]]$evidence$biological_model$species)       # Specie
              }
              else if (type == 'rna_expression'){
                if(length(response.disease$data[[k]]$literature$references[[1]]$lit_id>0)){
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$literature$references[[1]]$lit_id[1]) # paper url
                }else act.evidence$evidence <- c(act.evidence$evidence,'*')
                act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$comparison_name) # ej: "'oral squamous cell carcinoma' vs 'normal'"
              }
              else if (type == 'genetic_association'){
                if(length(response.disease$data[[k]]$literature$references[[1]]$list_id > 0))
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$literature$references[[1]]$lit_id) # paper url
                else act.evidence$evidence <- c(act.evidence$evidence,'*')
                
                if(length(response.disease$data[[k]]$evidence$gene2variant$functional_consequence)>0)
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$gene2variant$functional_consequence) # url ej: ontobee
                else act.evidence$evidence <- c(act.evidence$evidence,'*')
              }
              else if (type == 'somatic_mutation'){
                act.evidence$evidence <- c(act.evidence$evidence,
                                           response.disease$data[[k]]$evidence$known_mutations[[1]]$functional_consequence, # url ej: ontobee
                                           response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name) # ej: "sequence_alteration"
                if (length(response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer) > 0)
                  act.evidence$evidence <- c(act.evidence$evidence,response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer) # ej: "oncogene, TSG"
                else act.evidence$evidence <- c(act.evidence$evidence,'*')
              }
              else{
                act.evidence <- c()
                cat(paste('Unknown evidence type',response.disease$data[[k]]$type,'\n'))
              }
            }
            # Save actual evidence if is not empty
            if (length(act.evidence) > 0){
              if (! type %in% names(evidences)) evidences[[type]] <- c()
              evidences[[type]] <- c(evidences[[type]],list(act.evidence))
            }
          }
        }

        # Save found evidences for gene j (first row is empty)
        info[[geneList[j]]] <- evidences
      }else{
        # Not evidences found for gene j.
        # Add empty evidence
        info[[geneList[j]]] <- list()
      }
    }else{
      # Gene not found. 
      # Add empty evidence
      info[[geneList[j]]] <- list()
    }
  }
  
  cat("Evidences acquired successfully!\n")
  invisible(info)
}
