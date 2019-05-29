#' DEGsToDiseases obtains the information about what diseases are related to the DEGs indicated by parameter.
#'
#' The function obtains the information about what diseases are related to the DEGs indicated by parameter. For that, the function makes use of the web platform gene2Diseases.
#' @param geneList A list that contains the gene symbols or gene names of the DEGs.
#' @param minCitation Minimum number of citations of each genes in a disease to consider the genes related with the disease.
#' @return A list which contains the information about the diseases associated to each genes and to a sets of genes.
#' @examples
#' diseases <- DEGsToDiseases(c("KRT19","BRCA1"))

DEGsToDiseases <- function(geneList, minCitation = 5){

  if(length(geneList)[1] == 0 || is.null(geneList)){

    stop("The geneList is empty! Please, provide a right geneList.")

  }
  cat("Obtaining diseases related with the DEGs...\n")
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

  cat("Diseases acquired successfully!\n")
  results <- list(gene2Dis,result)
  names(results) <- c("DiseasesByGene","DiseasesByGeneSets")

  invisible(results)

}
