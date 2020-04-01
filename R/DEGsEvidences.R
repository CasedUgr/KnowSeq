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
  # Get disease id (it's necesary for evidences request)
  disease_ <- str_replace_all(disease,' ','-')
  r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease_,"&size=1&filter=disease",sep = ""))
  respon <- content(r_Ensembl)
  
  if ( 'size' %in% names(respon) && respon$size == 0){
    stop("Disease not found")
  }
  disease.id <- respon$data[[1]]$id

  if (verbose) cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
  base = "https://api.opentargets.io/v3/platform/public/evidence/filter?target="
  
  # Create empty output
  info <- list()
  
  # Iter in genes from geneList
  for(j in seq(length(unique(geneList)))){
    # Get gene id (it's necesary for evidences request)
    r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",geneList[j],"&size=1&filter=target",sep = ""))
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
          
          if ( grepl(disease,response.disease$data[[k]]$disease$efo_info$label)){
            # Create empty evidence
            act.evidence <- list()
            
            # Check if gene is matching
            if (response.disease$data[[k]]$target$gene_info$geneid == ensembl_id){
              type <- response.disease$data[[k]]$type
              # Save evidence score
              act.evidence <- list('score' = response.disease$data[[k]]$scores$association_score,
                                   'evidence' = list())
              # Save evidence codes
              act.evidence['code'] <- c()
              for (code.info in response.disease$data[[k]]$evidence$evidence_codes_info)
                act.evidence['codes'] <- c(code.info[[1]]$label)
              
              # Save information depending on evidence type
              if (type == 'known_drug'){
                act.evidence$evidence <- list(
                  'Drug Name'=response.disease$data[[k]]$drug$molecule_name, # Drug name
                  'Molecule Type'=response.disease$data[[k]]$drug$molecule_type) # Molecule type
              }
              else if (type == 'literature'){
                act.evidence$evidence <- list('Url'=response.disease$data[[k]]$evidence$literature_ref$lit_id) # paper url
              }
              else if(type == 'affected_pathway'){
                if (length(response.disease$data[[k]]$evidence$resource_score$method$reference) > 0){
                  act.evidence$evidence <- list('Url'=response.disease$data[[k]]$evidence$resource_score$method$reference) # paper url
                } else act.evidence$evidence <- list('Url'='*')
                if (length(response.disease$data[[k]]$evidence$known_mutations)>0){
                  act.evidence$evidence <- c(act.evidence$evidence,
                                             list('Name'=response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name))
                } else act.evidence$evidence <- c(act.evidence$evidence,list('Name'='*'))
                if(length( response.disease$data[[k]]$evidence$urls[[1]]$url ) > 0){
                  act.evidence$evidence <- c(act.evidence$evidence,
                                             list('Reactome Url'=response.disease$data[[k]]$evidence$urls[[1]]$url)) # reactome url
                } else act.evidence$evidence <- c(act.evidence$evidence,list('Reactome Url'='*'))
              }
              else if(type == 'animal_model'){
                act.evidence$evidence <- list(
                  'Is Associated'=response.disease$data[[k]]$evidence$biological_model$is_associated, # Boolean
                  'Specie'=response.disease$data[[k]]$evidence$biological_model$species)       # Specie
              }
              else if (type == 'rna_expression'){
                if(length(response.disease$data[[k]]$literature$references[[1]]$lit_id>0)){
                  act.evidence$evidence <- list('Url'=response.disease$data[[k]]$literature$references[[1]]$lit_id[1]) # paper url
                }else act.evidence$evidence <- list('Url'='*')
                act.evidence$evidence <- c(act.evidence$evidence,
                                           list('Comparison'=response.disease$data[[k]]$evidence$comparison_name)) # ej: "'oral squamous cell carcinoma' vs 'normal'"
              }
              else if (type == 'genetic_association'){
                if(length(response.disease$data[[k]]$literature$references[[1]]$list_id > 0))
                  act.evidence$evidence <- list('Url'=response.disease$data[[k]]$literature$references[[1]]$lit_id) # paper url
                else act.evidence$evidence <- list('Url'='*')
                
                if(length(response.disease$data[[k]]$evidence$gene2variant$functional_consequence)>0)
                  act.evidence$evidence <- c(act.evidence$evidence,
                                             list('Functional Consequence'=response.disease$data[[k]]$evidence$gene2variant$functional_consequence)) # url ej: ontobee
                else act.evidence$evidence <- c(act.evidence$evidence,list('Functional Consequence'='*'))
              }
              else if (type == 'somatic_mutation'){
                act.evidence$evidence <- list(
                  'Functional Consequence'=response.disease$data[[k]]$evidence$known_mutations[[1]]$functional_consequence, # url ej: ontobee
                  'Name'=response.disease$data[[k]]$evidence$known_mutations[[1]]$preferred_name) # ej: "sequence_alteration"
                if (length(response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer) > 0)
                  act.evidence$evidence <- c(act.evidence$evidence,
                                             list('Role in cancer'=response.disease$data[[k]]$evidence$known_mutations[[1]]$role_in_cancer)) # ej: "oncogene, TSG"
                else act.evidence$evidence <- c(act.evidence$evidence,list('Role in cancer'='*'))
              }
              else{
                act.evidence <- c()
                cat(paste('Unknown evidence type',response.disease$data[[k]]$type,'\n'))
              }
            }
            # Save actual evidence if it is not empty and if it is not repeated
            if (length(act.evidence) > 0){
              if (! type %in% names(evidences)){
                evidences[[type]] <- c()
                evidences[[type]] <- c(evidences[[type]],list(act.evidence))
              }
              else if (! any( unlist(lapply(list.map(evidences[[type]],evidence),function(x) identical(x,act.evidence$evidence)))))
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
    if (length(info[[geneList[j]]]) == 0) info[[geneList[j]]] = 'No  evidences found'
  }
  
  if (verbose) cat("Evidences acquired successfully!\n")
  invisible(info)
}