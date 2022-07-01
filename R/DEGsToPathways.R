#' The function uses the DEGs to retrieves the different pathways in which those DEGs involve any interaction.
#'
#' The function uses the DEGs to retrieves the different pathways in which those DEGs involve any interaction.
#' @param geneList A list which contains the DEGs that will be used to retrieve the related pathways to them.
#' @return A list with the pathways that contain relation to the DEGs within the geneList parameter.
#' @examples
#' DEGsToPathways(c("BRCA1","MLANA"))
#' 

DEGsToPathways <- function(geneList){
  
  dir <- system.file("extdata", package="KnowSeq")
  KEGG.data <- read.csv(paste(dir,"KEGGPathsDB.csv",sep = "/"))
  Gene.data <- read.csv(paste(dir,"Genes_length_Homo_Sapiens.csv",sep = "/"))
  KEGG.matrix <- matrix(ncol = 5)
  cat("Retrieving information about KEGG pathways...\n")
  
  for(gene in geneList){
    if(!is.na(gene)){
      pos <- which(gene == Gene.data$Gene_name)
      pathways <- Gene.data$keggPathways[pos]
      pathways <- unlist(strsplit(pathways, split = ", "))
      
      for(pathway in pathways){
        if (length(pathway) == 0) break
          path_map <- pathway
          if(length(which(path_map == KEGG.matrix[,1])) == 0){
            KEGG.matrix <- rbind(KEGG.matrix,c(KEGG.data[which(path_map == KEGG.data$KEGG_Id),],gene))
          }else if(length(which(path_map == KEGG.matrix[,1])) > 0){
            pos <- which(path_map == KEGG.matrix[,1])
            KEGG.matrix[pos,5] <- paste(KEGG.matrix[pos,5], gene, sep = ", ")
          }
      }
      
    }
  }
  
  KEGG.matrix <- as.data.frame(KEGG.matrix,stringsAsFactors=FALSE)
  KEGG.matrix <- KEGG.matrix[-1,]
  names(KEGG.matrix) <- c("KEGG_Path","Name","Description","Class","Genes")

  rownames(KEGG.matrix) <- NULL
  
  invisible(KEGG.matrix)

}
