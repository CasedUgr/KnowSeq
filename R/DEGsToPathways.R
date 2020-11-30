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
  
  KEGG.data <- matrix(ncol = 5)
  
  cat("Retrieving information about KEGG pathways...\n")
  
  for(gene.index in seq(dim(myAnnotation)[1])){
    gene <- myAnnotation[gene.index,'entrezgene_id']
    gene.name <- myAnnotation[gene.index,'external_gene_name']
    if(!is.na(gene)){
      get_GO <- GET(paste("http://rest.kegg.jp/get/hsa:",gene,sep = ""))
      get_GO_text <- content(get_GO, "text")
      pathway_start <- str_locate_all(pattern = "PATHWAY", get_GO_text)[[1]][2]
      pathway_end <- str_locate_all(pattern = "BRITE", get_GO_text)[[1]][1]
      pathways <- substr(get_GO_text,pathway_start+1,pathway_end-1)
      pathways <- strsplit(pathways,split = "\n")
      
      index <- grep("hsa", unlist(pathways))
      
      for(i in index){
        pathway <- as.character(unlist(str_extract_all(unlist(pathways)[i],"hsa[a-zA-Z0-9]{5}")))
        if (length(pathway) == 0) break
        get_pathway <- GET(paste("http://rest.kegg.jp/get/map",substr(pathway,4,nchar(pathway)),sep = ""))
        get_pathway_text <- content(get_pathway, "text")
      
        start <- str_locate_all(pattern = "ENTRY", get_pathway_text)[[1]][2]
        final <- str_locate_all(pattern = "Pathway", get_pathway_text)[[1]][2]
        path.index <- substr(get_pathway_text,start,final)
        path.index <- strsplit(path.index,split = " ")[[1]]
        path <- path.index[path.index != ""]
        path <- path[2]
        
        start <- str_locate_all(pattern = "\nNAME ", get_pathway_text)[[1]][2]
        final <- str_locate_all(pattern = "\nDESCRIPTION", get_pathway_text)[[1]][2]
        name <- substr(get_pathway_text,start + 8,final - 12)
        
        start <- str_locate_all(pattern = "\nDESCRIPTION ", get_pathway_text)[[1]][2]
        final <- str_locate_all(pattern = "\nCLASS", get_pathway_text)[[1]][2]
        description <- substr(get_pathway_text,start + 1,final - 6)
        description <- gsub("\"","",description)
        
        start <- str_locate_all(pattern = "\nCLASS ", get_pathway_text)[[1]][2]
        final <- str_locate_all(pattern = "\nPATHWAY_MAP", get_pathway_text)[[1]][2]
        class <- substr(get_pathway_text,start + 7,final - 12)

        KEGG.data <- rbind(KEGG.data,c(as.character(path),as.character(name), as.character(description),as.character(class),as.character(gene.name)))
      }
      
    }
  }
  KEGG.data <- KEGG.data[-1,]
  KEGG.data <- as.data.frame(KEGG.data,stringsAsFactors=FALSE)
  names(KEGG.data) <- c("KEGG_Path","Name","Description","Class","Genes")
  naPos <- which(is.na(KEGG.data) == TRUE)
  
  # Collapse repeated pathways
  remove.index <- c()
  repeated.pathways <- names(which(table(KEGG.data$KEGG_hsa) > 1))
  for (pathway in repeated.pathways){
    index <- which(KEGG.data$KEGG_hsa == pathway,)
    genes <- KEGG.data[index,'Genes']
    KEGG.data[index[1],'Genes'] <- paste(genes,collapse=', ')
    remove.index <- c(remove.index,index[-1])
  }
  
  if(!is.null(remove.index))
    KEGG.data <- KEGG.data[-remove.index,]
  
  rownames(KEGG.data) <- NULL
  
  
  #cat("Retrieving information about Reactome pathways...\n")
  
  response <- GET("https://reactome.org/content/query?q=BRCA1&species=Homo+sapiens&species=Entries+without+species&cluster=true")
  response_text <- content(response, "text")
  proteinMatches <- str_locate_all(pattern = "Protein", response_text)
  proteinEnd <- unlist(proteinMatches)[length(unlist(proteinMatches))]
  geneMatches <- str_locate_all(pattern = "BRCA1", response_text)
  geneMatches <- unlist(geneMatches)
  geneInfo <- substr(response_text,proteinEnd, geneMatches[which(geneMatches > proteinEnd)[1]])
  geneInfo <- regmatches(geneInfo,gregexpr("(R-HSA-\\d{1,})", geneInfo, perl=TRUE))
  
  invisible(KEGG.data)

}
