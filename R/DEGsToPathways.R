#' The function uses the DEGs to retrieves the different pathways in which those DEGs involve any interaction.
#'
#' The function uses the DEGs to retrieves the different pathways in which those DEGs involve any interaction.
#' @param geneList A list which contains the DEGs that will be used to retrieve the related pathways to them.
#' @return A list with the pathways that contain relation to the DEGs within the geneList parameter.
#' @examples
#' DEGsToPathways(c("BRCA1","MLANA"))
#' 

DEGsToPathways <- function(geneList){
  
  myAnnotation <- getGenesAnnotation(geneList,attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id"),filter="external_gene_name",notHSapiens = FALSE)
  
  myAnnotation <- myAnnotation[which(!is.na(myAnnotation$entrezgene_id) == TRUE),]
  dir <- system.file("extdata", package="KnowSeq")
  KEGG.data <- read.csv(paste(dir,"KEGGPathsDB.csv",sep = "/"))
  KEGG.matrix <- matrix(ncol = 5)
  cat("Retrieving information about KEGG pathways...\n")
  
  for(gene.index in seq(dim(myAnnotation)[1])){
    gene <- myAnnotation[gene.index,'entrezgene_id']
    gene.name <- myAnnotation[gene.index,'external_gene_name']
    if(!is.na(gene)){
      get_GO <- GET(paste("https://rest.kegg.jp/get/hsa:",gene,sep = ""))
      get_GO_text <- content(get_GO, "text")
      pathway_start <- str_locate_all(pattern = "PATHWAY", get_GO_text)[[1]][2]
      pathway_end <- str_locate_all(pattern = "BRITE", get_GO_text)[[1]][1]
      pathways <- substr(get_GO_text,pathway_start+1,pathway_end-1)
      pathways <- strsplit(pathways,split = "\n")
      
      index <- grep("hsa", unlist(pathways))
      
      for(i in index){
        pathway <- as.character(unlist(str_extract_all(unlist(pathways)[i],"hsa[a-zA-Z0-9]{5}")))
        if (length(pathway) == 0) break
          path_map <- gsub("hsa","map",pathway)
          if(length(which(path_map == KEGG.matrix[,1])) == 0){
            KEGG.matrix <- rbind(KEGG.matrix,c(KEGG.data[which(path_map == KEGG.data$KEGG_Id),],gene.name))
          }else if(length(which(path_map == KEGG.matrix[,1])) > 0){
            pos <- which(path_map == KEGG.matrix[,1])
            KEGG.matrix[pos,5] <- paste(KEGG.matrix[pos,5], gene.name, sep = ", ")
          }
      }
      
    }
  }
  
  KEGG.matrix <- as.data.frame(KEGG.matrix,stringsAsFactors=FALSE)
  KEGG.matrix <- KEGG.matrix[-1,]
  names(KEGG.matrix) <- c("KEGG_Path","Name","Description","Class","Genes")

  rownames(KEGG.matrix) <- NULL
  
  
  # For future implementation of Reactome
  #cat("Retrieving information about Reactome pathways...\n")
  
  #response <- GET("https://reactome.org/content/query?q=BRCA1&species=Homo+sapiens&species=Entries+without+species&cluster=true")
  #response_text <- content(response, "text")
  #proteinMatches <- str_locate_all(pattern = "Protein", response_text)
  #proteinEnd <- unlist(proteinMatches)[length(unlist(proteinMatches))]
  #geneMatches <- str_locate_all(pattern = "BRCA1", response_text)
  #geneMatches <- unlist(geneMatches)
  #geneInfo <- substr(response_text,proteinEnd, geneMatches[which(geneMatches > proteinEnd)[1]])
  #geneInfo <- regmatches(geneInfo,gregexpr("(R-HSA-\\d{1,})", geneInfo, perl=TRUE))
  
  invisible(KEGG.matrix)

}
