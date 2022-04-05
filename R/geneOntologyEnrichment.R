#' geneOntologyEnrichment obtains the information about what Gene Ontology terms are related to the DEGs.
#'
#' The function obtains the information about GO terms from the three differents ontologies that are related to the DEGs. The function also returns the description about each GO and a list of genes that are inside of each GO.
#' @param geneList A list that contains entrez gene id of the DEGs. Entrez gene id can be obtained using getAnnotationFromEnsembl function.
#' @param ontologies A list that contains ontologies to be searchs. Values must be contained in the following three: BP, CC, MF.
#' @param pvalCutOff The maximum p-value to considers that a genes is related with a GO term.
#' @param geneType A string indicating the type of genes in geneList, it must be one of indicated in DAVIDs API documentation.
#' @return A list that contains a matrix for each of the possible ontologies and a matrix with the GOs for the three ontologies together.
#' @examples
#' \dontrun{GOsList <- geneOntologyEnrichment(data$entrezgene_id,geneType='ENTREZ_GENE_ID',pvalCutOff=0.1)}

geneOntologyEnrichment <- function(geneList, geneType="ENTREZ_GENE_ID", ontologies=c('BP','CC','MF'), pvalCutOff=1){
  if(!is(geneList)[1]=='character'){stop('The class of geneList must be character')}
  if(!is(geneType)[1]=='character'){stop('The class of geneType must be character')}
  if(!is(ontologies)[1]=='character'){stop('The class of ontologies must be character')}
  if(!is.numeric(pvalCutOff)){stop("The class of pvalCutOff parameter must be numeric.")}

  annotations=''
  for ( ontology in ontologies){
    if (ontology=='BP') annotations <- paste(annotations,',','GOTERM_BP_ALL',sep='')
    else if (ontology=='CC') annotations <- paste(annotations,',','GOTERM_CC_ALL',sep='')
    else if (ontology=='MF') annotations <- paste(annotations,',','GOTERM_MF_ALL',sep='')
    else stop(paste('Ontology',ontology,'not found. Ontology values must be contained in the following three: BP, CC, MF'))
  }
  annotations <- substr(annotations,2,nchar(annotations))
  
  if ( geneType == 'ENTREZ_GENE_ID')
    {gene.type <- 'entrezgene_id'}
  else if (geneType == 'GENE_SYMBOL')
    {gene.type <- 'external_gene_name'}
  else 
    {gene.type <- tolower(geneType)}
  
  cat('Getting gene symbols...')
  genes.annotations <- getGenesAnnotation(geneList,attributes=c("external_gene_name","entrezgene_id"),filter=gene.type)
  

  cat('Retrieving Gene Ontology terms related to the list of DEGs...\n')
  geneList <- paste(genes.annotations$entrezgene_id, collapse=",")
  base  <- 'https://david.ncifcrf.gov/'
  
  url <- paste(base,'api.jsp?type=','ENTREZ_GENE_ID','&ids=',geneList ,'&tool=chartReport&annot=',annotations, sep='') # Do not change order
  
  response <- GET(url)
  response <- content(response,'text', encoding = "utf-8")
  
  rowids <- str_match(response,'document.apiForm.rowids.value=\"(.*?)"')[2]
  annotids <- str_match(response,'document.apiForm.annot.value=\"(.*?)"')[2]

  url <- paste(base,'chartReport.jsp?rowids=', rowids, '&annot=', annotids,'&count=0&ease=',pvalCutOff, sep='')
  response <- GET(url)
  response <- content(response,'text', encoding = "utf-8")
  
  downloadFileName <- str_match(response,'<a href="data/download/(.*?)"')[2]
  
  if (grepl(".txt", downloadFileName)){
    downloadFileName <- paste(base,"data/download/", downloadFileName, sep="")
  }else{
    stop('Error in request. Please check parameters or use geneType as ENTREZ_GENE_ID or GENE_SYMBOL')
  }

  response <- GET(downloadFileName)
  response <- content(response)
  
  index <- gregexpr(pattern ='Category',response)[[1]][1]
  response <- substr(response,index,nchar(response))
  gos.data <- read.table(text=response,sep='\t',header=TRUE)

  final.gos <- list()
  for (go.type in c('GOTERM_MF_ALL','GOTERM_CC_ALL','GOTERM_BP_ALL')){
    act.gos <- gos.data[gos.data$Category == go.type,colnames(gos.data)!='Category']
    if(dim(act.gos)[1]>0){
      ontology.term <- c()
      remove <- c()
      for ( i in seq(length(act.gos$Term))){
        ontology.str <- as.character(act.gos[i,'Term'])
        num.gos <- str_count(ontology.str,'GO:')
        if (num.gos == 1) ontology.term <- c(ontology.term,strsplit(ontology.str,'~'))
        else remove <- c(remove,i)
      }
      ontology.term <- t(matrix(unlist(ontology.term),nrow=2))
      
      if (length(remove)>0) act.gos <- act.gos[-remove,]
      act.gos[['GO.ID']] <- ontology.term[,1]
      act.gos[['Term']] <- ontology.term[,2] 

      # Remove NAs
      genes.annotations <- genes.annotations[!is.na(genes.annotations$entrezgene_id),]
      
      # Add column with gene symbols
      gene.names <- act.gos$Genes
      for ( gen in seq(dim(genes.annotations)[1])){
        gene.names <- str_replace(gene.names,as.character(genes.annotations[gen,'entrezgene_id']),
                                  as.character(genes.annotations[gen,'external_gene_name']))
      }
      act.gos['Gene Symbols'] <- gene.names
      

      descriptions <- c()
      for (go in act.gos[['GO.ID']]){
        go <- gsub(':','%3A',go)
        get_GO <- GET(paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",go,sep=''))
        get_GO_text <- content(get_GO,'text', encoding = "utf-8")
        get_GO_json <- fromJSON(get_GO_text, flatten = TRUE)
        descriptions  <- c(descriptions,as.character(get_GO_json$results$definition.text))
      }
      act.gos[['Description']] <- descriptions
      act.gos <- act.gos[,c("GO.ID",colnames(act.gos)[colnames(act.gos)!='GO.ID'])]
    }
    else{
      print(paste('Empty',go.type))
    }
    index <- gregexpr(pattern ='_',go.type)[[1]]
    onto.name <- substr(go.type,index[1]+1,index[2]-1)
    final.gos[[paste(onto.name,'Ontology GOs')]] <- act.gos
  }
  final.gos[['All Ontologies GO']] <- rbind(final.gos[['MF Ontology GOs']],final.gos[['CC Ontology GOs']],final.gos[['BP Ontology GOs']])
  
  return(final.gos)
}

