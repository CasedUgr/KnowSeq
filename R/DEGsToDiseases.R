#' DEGsToDiseases obtains the information about what diseases are related to the DEGs indicated by parameter.
#'
#' The function obtains the information about what diseases are related to the DEGs indicated by parameter. For that, the function makes use of the web platforms gene2Diseases and targetValidation.
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param minCitation Minimum number of citations of each genes in a disease to consider the genes related with the disease.
#' @param size The number of diseases to retrieve from targetValidation
#' @param method The name of the desired web platform to use for the diseases download: genes2Diseases or targetValidation
#' @return A list which contains the information about the diseases associated to each genes or to a set of genes.
#' @examples
#' diseases <- DEGsToDiseases(c("KRT19","BRCA1"))

DEGsToDiseases <- function(geneList, minCitation = 5, size = 10, method = "targetValidation"){

  if(length(geneList)[1] == 0 || is.null(geneList)){

    stop("The geneList is empty! Please, provide a right geneList.")

  }
  
  if(method == "genes2Diseases"){
    
      cat("Obtaining related diseases with the DEGs from genes2Diseases platform...\n")
      base = "http://cbdm-01.zdv.uni-mainz.de/~jfontain/cgi-bin/genes2diseases.pl"
    
      genePeti = character()
    
      for(i in seq_len(length(geneList))){
    
        genePeti <- paste(genePeti,geneList[i],sep = "|")
    
      }
    
      genePeti <- substr(genePeti,2,nchar(genePeti))
    
      petition = paste("?analysis_type=gene&items=",genePeti,"&min_citations_for_gene=", minCitation,"&fdr_cutoff=0.1&output_type=text",sep = "")
    
      download.file(url = paste(base,petition,sep = ""),destfile="gene2diseases.tsv",method="libcurl",quiet = TRUE)
      gene2Dis <- read.csv("gene2diseases.tsv",sep = "\t")
    
      file.remove("gene2diseases.tsv")
    
      geneSetPeti = character()
    
      for(i in seq_len(length(geneList))){
    
        geneSetPeti <- paste(geneSetPeti,geneList[i],sep = "|")
    
      }
    
      geneSetPeti <- substr(geneSetPeti,2,nchar(geneSetPeti))
    
      petition = paste("?analysis_type=geneset&items=",geneSetPeti,"&min_citations_for_gene=", minCitation,"&fdr_cutoff=0.1&output_type=text",sep = "")
    
      download.file(url = paste(base,petition,sep = ""),destfile="geneset2diseases.tsv",method="libcurl",quiet = TRUE)
      geneset2Dis <- read.csv("geneset2diseases.tsv",sep = "\t")
    
      file.remove("geneset2diseases.tsv")
    
    
      cat("Obtaining diseases related with sets of DEGs...\n")
    
      result <- vector("list", length(unique(geneset2Dis$Gene.symbols)))
      names(result) <- unique(geneset2Dis$Gene.symbols)
    
      for(i in seq_len(length(unique(geneset2Dis$Gene.symbols)))){
    
        currentGene <- geneset2Dis[which(unique(geneset2Dis$Gene.symbols)[i] == geneset2Dis$Gene.symbols),]
        currentGene$Genes.count <- NULL
        currentGene$Genes.percentage <- NULL
        currentGene$Gene.symbols <- NULL
        currentGene$Entrez.Gene.IDs <- NULL
        result[[i]] <- currentGene
    
      }
    
      results <- list(gene2Dis,result)
      names(results) <- c("DiseasesByGene","DiseasesByGeneSets")
      
  }else if(method == "targetValidation"){
    
      cat("Obtaining related diseases with the DEGs from targetValidation platform...\n")
      base = "https://api.opentargets.io/v3/platform/public/association/filter?target="
      
      genePeti = character()
      results <- vector("list", 0)
      
      for(j in seq_len(length(unique(geneList)))){
        r_Ensembl <- httr::GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",geneList[j],"&size=1&filter=target",sep = ""))
        respon <- content(r_Ensembl)
        ensembl_id <- respon$data[[1]]$data$ensembl_gene_id
        
        if(!is.na(ensembl_id)){
          
            genePeti = paste(base, ensembl_id,"&direct=true&fields=association_score&fields=disease.efo_info.label&size=",size,sep = "")
            r <- GET(genePeti)
            response = content(r)
            
            if(response$size != 0){
                currentGene = matrix(, nrow = response$size, ncol = 9)
                
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
                  
                }
                colnames(currentGene) <- c("Disease","Overall Score","Literature","RNA Expr.","Genetic Assoc.","Somatic Mut.","Known Drug","Animal Model","Affected Pathways")
                results[[j]] <- currentGene
                names(results)[j] <- unique(geneList)[j]
            }else{
              
              cat(paste("Removing ", unique(geneList)[j], " from the list, there is no exists associated diseases for this gene.\n",sep = ""))
              
            }
          }
        }
  
  }
  
  cat("Diseases acquired successfully!\n")
  results <- results[lengths(results) != 0]
  invisible(results)

}
