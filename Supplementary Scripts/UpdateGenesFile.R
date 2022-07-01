require(httr)
require(stringr)
require(KnowSeq)

dir <- system.file("extdata", package="KnowSeq")
Gene.data <- read.csv(paste(dir,"Genes_length_Homo_Sapiens.csv",sep = "/"))
Gene.info <- read.csv(paste(dir,"GRCh38Annotation.csv",sep = "/"))
matches <- match(Gene.data$Gene_name,Gene.info$external_gene_name)

keggPathways <- c()

for(gene.index in Gene.data$Gene_name){
  if(length(which(gene.index == Gene.info$external_gene_name)) > 0){
    pos <- which(gene.index == Gene.info$external_gene_name)
    gene <- Gene.info$entrezgene_id[pos]
    gene.name <- gene.index
    gene = gene[1]
    path_gene = ""
    
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
        if (length(pathway) != 0){
          path_map <- gsub("hsa","map",pathway)
          if(path_gene == "")
            path_gene <- path_map
          else
            path_gene <- paste(path_gene, path_map, sep = ", ")
        }
      }
    }
  }
  keggPathways <- c(keggPathways,path_gene)
}


Gene.data <- cbind(Gene.data,keggPathways)
write.csv(Gene.data,"inst/extdata/Genes_length_Homo_Sapiens.csv", row.names = F)
